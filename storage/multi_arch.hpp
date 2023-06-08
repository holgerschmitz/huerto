/*
 * multi_arch.hpp
 *
 * Created on: 08 Jun 2023
 * Author: Holger Schmitz (holger@notjustphysics.com)
 */

#include <exception>
#include <type_traits>

/**
 * @brief Store data across multiple architectures and access it from host and device
 * 
 * @tparam ArchitecturePolicy The policy that determines where data is stored
 */
template<class ArchitecturePolicy>
class MultiArchData {
  public:
    /**
     * @brief The type of data stored on the host 
     */
    typedef ArchitecturePolicy::HostDataType HostDataType;

    /**
     * @brief The type of data stored on the device 
     */
    typedef ArchitecturePolicy::DeviceDataType DeviceDataType;
  private:

    /**
     * @brief The data stored on the host
     */
    HostDataType hostData;

    /**
     * @brief The data stored on the device
     */
    DeviceDataType deviceData;

    /**
     * @brief An enum indicating where the data was last accessed
     */
    enum { host, device, unset } lastAccess = unset;
  public:
    /**
     * @brief Set the host data
     * 
     * @param hostData the data to be stored on the host
     */
    void setHostData(HostDataType hostData);

    /**
     * @brief Access the host data
     * 
     * This function will copy the data to the host if it was last accessed from the device
     * 
     * @return HostDataType The data stored on the host
     */
    HostDataType accessHostData();

    /**
     * @brief Access the device data
     * 
     * This function will copy the data to the device if it was last accessed from the host
     * 
     * @return DeviceDataType The data stored on the device
     */
    DeviceDataType accessDeviceData();
};

/**
 * @brief A policy that defines host and device data to be the same
 * 
 * Copy operations are no-ops and simply return the input
 * 
 * @tparam DataType The data type to be stored
 */
template<class DataType>
struct SingleArchDataPolicy {
    /**
     * @brief The type of data stored on the host 
     */
    typedef DataType HostDataType;
    
    /**
     * @brief The type of data stored on the device 
     */
    typedef DataType DeviceDataType;

    /**
     * @brief No-op copy to device
     * 
     * The input is simply returned
     * 
     * @param hostData The data to be copied to the device
     * @return DeviceDataType The data stored on the device
     */
    static DeviceDataType copyToDevice(HostDataType hostData);

    /**
     * @brief No-op copy to host
     * 
     * The input is simply returned
     * 
     * @param deviceData The data to be copied to the host
     * @return HostDataType The data stored on the host
     */
    static HostDataType copyToHost(DeviceDataType deviceData);
};

//=================================================================
//=============== SingleArchDataPolicy ============================
//=================================================================


template<class DataType>
SingleArchDataPolicy<DataType>::DeviceDataType SingleArchDataPolicy<DataType>::copyToDevice(HostDataType hostData)
{
  return hostData;
}

template<class DataType>
SingleArchDataPolicy<DataType>::HostDataType SingleArchDataPolicy<DataType>::copyToHost(DeviceDataType deviceData)
{
  return deviceData;
}

template<class ArchitecturePolicy>
void MultiArchData<ArchitecturePolicy>::setHostData(HostDataType hostData)
{
  this->hostData = hostData;
}

template<class ArchitecturePolicy>
MultiArchData<ArchitecturePolicy>::HostDataType MultiArchData<ArchitecturePolicy>::accessHostData()
{
  if (lastAccess == unset) {
    throw std::exception("MultiArchData no initialised before access");
  }
  if (lastAccess == device) {
    hostData = ArchitecturePolicy::copyToHost(deviceData);
  }
  lastAccess = host;
  return hostData;
}

template<class ArchitecturePolicy>
MultiArchData<ArchitecturePolicy>::DeviceDataType MultiArchData<ArchitecturePolicy>::accessDeviceData()
{
  if (lastAccess == unset) {
    throw std::exception("MultiArchData no initialised before access");
  }
  if (lastAccess == host) {
    deviceData = ArchitecturePolicy::copyToDevice(hostData);
  }
  lastAccess = device;
  return deviceData;
}