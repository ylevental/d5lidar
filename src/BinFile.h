#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

class BinFile {
 public:
  BinFile() = default;

  BinFile(const std::string& filename) { open(filename); }

  ~BinFile() { close(); }

  void open(const std::string& filename);

  void close();

  struct FileHeader {
    /// The "magic string" used to verify the file type ("DIRSIGPROTO")
    char fileIdentifier[11];

    /// The version of file format specification.
    char fileFormatRev;

    /// The byte order, 0 indicates big endian, 1 for little endian
    char byteOrder;

    /// The date/time that file was created. Format: YYYYMMDDhhmm.ss
    char creationDateTime[15];

    /// The DIRSIG version used to create file.
    char dirsigVersion[32];

    /// The arbitrary user-supplied text describing simulation.
    char simulationDesc[256];

    /// The latitude of scene origin in degrees using +East [-180,180)
    double sceneOriginLatitude = 0;

    /// The longitude of scene origin in degrees using +North [-90,90]
    double sceneOriginLongitude = 0;

    /// The height of scene origin above WGS84 reference ellipsoid [m]
    double sceneOriginHeight = 0;

    /// The description of mount type for laser
    char laserMountType[16];

    /// The description of mount type for detector
    char detectorMountType[16];

    /// The number of detector elements in the X dimension
    uint32_t xDetectorCount = std::numeric_limits<uint32_t>::max();

    /// The number of detector elements in the Y dimension
    uint32_t yDetectorCount = std::numeric_limits<uint32_t>::max();

    /// The distance between detector centers [microns] in X dimension
    double xDetectorPitch = -1;

    /// The distance between detector centers [microns] in Y dimension
    double yDetectorPitch = -1;

    /// The distance between array center and optical axis [microns]
    double xArrayOffset = 0;

    /// The distance between array center and optical axis [microns]
    double yArrayOffset = 0;

    /// The first radial lens distortion coefficient k1
    double lensK1 = 0;

    /// The second radial lens distortion coefficient k2
    double lensK2 = 0;

    /// The number of tasks in file
    uint32_t taskCount = 0;

    /// The ID for the focal plane array
    uint16_t fpaId = 0;
  };

  struct TaskHeader {
    /// The ASCII/Text description of the task
    char description[64];

    /// The date/time for task start. Format: YYYYMMDDhhmm.ss
    char startDateTime[15];

    /// The date/time for task stop. Format: YYYYMMDDhhmm.ss
    char stopDateTime[15];

    /// The instrument focal length in [mm]
    double focalLength = 0;

    /// The nominal laser pulse repetition frequency [Hz]
    double pulseRepetitionFreq = 0;

    /// The nominal laser pulse duration [sec]
    double pulseDuration = 0;

    /// The nominal laser pulse energy [J]
    double pulseEnergy = 0;

    /// The laser spectral center[microns]
    double laserSpectralCenter = 0;

    /// The laser spectral width [microns]
    double laserSpectralWidth = 0;

    /// The number of pulses in this task
    uint32_t pulseCount = 0;
  };

  struct PulseHeader {
    /// The pulse time relative to task start [sec]
    double pulseTime = -1;

    /// The time gating start time relative to task start [sec]
    double timeGateStart = -1;

    /// The time gating stop time relative to task start [sec]
    double timeGateStop = -1;

    /// The number of temporal samples (equally spaced) over time gating
    uint32_t timeGateBinCount = 0;

    /// The number of samples (equally spaced) within a time bin
    uint32_t samplesPerTimeBin = 1;

    /// The platform location (x,y,z) in scene local coordinates [m]
    double platformLocation[3];

    /// The platform Euler angles (x,y,z) applied in "yzx" order [rad]
    double platformRotation[3];

    /// The Tx to mount transform (4x4 affine matrix, row major)
    double transmitterToMount[16];

    /// The Tx pointing Euler angles (x,y,z) applied in "yzx" order [rad]
    double transmitterMountPointing[3];

    /// The Tx mount to platform transform (4x4 matrix, row major)
    double transmitterMountToPlatform[16];

    /// The Rx to mount transform (4x4 affine matrix, row major)
    double receiverToMount[16];

    /// The Rx pointing Euler angles (x,y,z) applied in "yzx" order [rad]
    double receiverMountPointing[3];

    /// The Rx mount to platform transform (4x4 matrix, row major)
    double receiverMountToPlatform[16];

    /**
     * Data type used in pulse results with only type 5 (double precision)
     * being currently supported. Corresponds to the ENVI type numbers.
     */
    int pulseDataType = 5;

    /**
     * The compression mode.
     * 0 indicates no compression,
     * 1 indicates zlib compression of pulse data
     */
    char dataCompression = 0;

    /// The pulse index (starts at 0) within a task
    uint32_t pulseIndex = 0;

    /// The number of bytes of (possibly compressed) pulse data
    uint64_t pulseDataBytes = 0;

    /// The system polarization, transmit (optics, scan mirrors, etc...)
    double systemPolTransmit[16];

    /// The system polarization, receive (optics, scan mirrors, etc...)
    double systemPolReceive[16];
  };

  FileHeader fileHeader;

  struct Task : TaskHeader {
    std::vector<size_t> pulseOffsets;
  };

  std::vector<Task> tasks;

  PulseHeader loadPulseHeader(const Task& task, size_t pulseIndex);

  std::vector<double> loadPulse(const Task& task, const PulseHeader& pulse);

 private:
  bool swapEndian_ = false;

  std::ifstream fileStream_;

  void fileRead(void* ptr, size_t num, size_t sz);

  void fileRead(TaskHeader& header);

  void fileRead(PulseHeader& header);

  template <
      typename T,
      typename = std::enable_if_t<std::is_arithmetic_v<T>, void>>
  void fileRead(T& value) {
    fileRead(&value, 1, sizeof(T));
  }

  template <typename T, size_t N>
  void fileRead(T (&value)[N]) {
    fileRead(&value[0], N, sizeof(T));
  }

  template <typename T, typename... U>
  void fileReadMany(T&& first, U&&... others) {
    (fileRead(std::forward<T>(first)));
    (fileRead(std::forward<U>(others)), ...);
  }
};
