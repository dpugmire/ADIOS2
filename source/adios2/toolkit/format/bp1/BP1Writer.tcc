/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP1Writer.tcc
 *
 *  Created on: Apr 11, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */
#ifndef ADIOS2_TOOLKIT_FORMAT_BP1_BP1WRITER_TCC_
#define ADIOS2_TOOLKIT_FORMAT_BP1_BP1WRITER_TCC_

#include "BP1Writer.h"

#include "adios2/helper/adiosFunctions.h"

namespace adios
{
namespace format
{

// PUBLIC
template <class T>
BP1Writer::ResizeResult BP1Writer::ResizeBuffer(const Variable<T> &variable)
{
    size_t variableData =
        GetVariableIndexSize(variable) + variable.PayLoadSize();
    size_t requiredCapacity = variableData + m_HeapBuffer.m_DataPosition;

    if (requiredCapacity > m_MaxBufferSize && m_MaxBufferSize > 0) // is set
    {
        if (m_HeapBuffer.GetDataSize() < m_MaxBufferSize)
        {
            m_HeapBuffer.ResizeData(m_MaxBufferSize);
            return ResizeResult::FLUSH;
        }
    }

    return ResizeResult::SUCCESS;
}

template <class T>
void BP1Writer::WriteVariableMetadata(const Variable<T> &variable) noexcept
{
    if (m_Profiler.IsActive)
    {
        m_Profiler.Timers.at("buffering").Resume();
    }

    Stats<typename TypeInfo<T>::ValueType> stats = GetStats(variable);

    stats.TimeIndex = m_MetadataSet.TimeStep;
    // Get new Index or point to existing index
    bool isNew = true; // flag to check if variable is new
    BP1Index &variableIndex =
        GetBP1Index(variable.m_Name, m_MetadataSet.VarsIndices, isNew);
    stats.MemberID = variableIndex.MemberID;

    // write metadata header in data and extract offsets
    stats.Offset = m_HeapBuffer.m_DataAbsolutePosition;
    WriteVariableMetadataInData(variable, stats);
    stats.PayloadOffset = m_HeapBuffer.m_DataAbsolutePosition;

    // write to metadata  index
    WriteVariableMetadataInIndex(variable, stats, isNew, variableIndex);

    ++m_MetadataSet.DataPGVarsCount;

    if (m_Profiler.IsActive)
    {
        m_Profiler.Timers.at("buffering").Pause();
    }
}

template <class T>
void BP1Writer::WriteVariablePayload(const Variable<T> &variable) noexcept
{
    if (m_Profiler.IsActive)
    {
        m_Profiler.Timers.at("buffering").Resume();
    }

    CopyToBufferThreads(m_HeapBuffer.m_Data, m_HeapBuffer.m_DataPosition,
                        variable.m_AppValues, variable.TotalSize(), m_Threads);
    m_HeapBuffer.m_DataAbsolutePosition += variable.PayLoadSize();

    if (m_Profiler.IsActive)
    {
        m_Profiler.Timers.at("buffering").Pause();
    }
}

// PRIVATE FUNCTIONS
template <class T>
size_t BP1Writer::GetVariableIndexSize(const Variable<T> &variable) const
    noexcept
{
    // size_t indexSize = varEntryLength + memberID + lengthGroupName +
    // groupName + lengthVariableName + lengthOfPath + path + datatype
    size_t indexSize = 23; // without characteristics
    indexSize += variable.m_Name.size();

    // characteristics 3 and 4, check variable number of dimensions
    const size_t dimensions = variable.m_Count.size();
    indexSize += 28 * dimensions; // 28 bytes per dimension
    indexSize += 1;               // id

    // characteristics, offset + payload offset in data
    indexSize += 2 * (1 + 8);
    // characteristic 0, if scalar add value, for now only allowing string
    if (dimensions == 1)
    {
        indexSize += sizeof(T);
        indexSize += 1; // id
        // must have an if here
        indexSize += 2 + variable.m_Name.size();
        indexSize += 1; // id
    }

    // characteristic statistics
    if (m_Verbosity == 0) // default, only min and max
    {
        indexSize += 2 * (sizeof(T) + 1);
        indexSize += 1 + 1; // id
    }

    return indexSize + 12; // extra 12 bytes in case of attributes
    // need to add transform characteristics
}

template <class T>
BP1Writer::Stats<typename TypeInfo<T>::ValueType>
BP1Writer::GetStats(const Variable<T> &variable) const noexcept
{
    Stats<typename TypeInfo<T>::ValueType> stats;
    const std::size_t valuesSize = variable.TotalSize();

    if (m_Verbosity == 0)
    {
        GetMinMaxThreads(variable.m_AppValues, valuesSize, stats.Min, stats.Max,
                         m_Threads);
    }
    return stats;
}

template <class T>
void BP1Writer::WriteVariableMetadataInData(
    const Variable<T> &variable,
    const Stats<typename TypeInfo<T>::ValueType> &stats) noexcept
{
    auto &buffer = m_HeapBuffer.m_Data;
    auto &position = m_HeapBuffer.m_DataPosition;

    // for writing length at the end
    const size_t varLengthPosition = position;
    position += 8; // skip var length (8)

    CopyToBuffer(buffer, position, &stats.MemberID);

    WriteNameRecord(variable.m_Name, buffer, position);
    position += 2; // skip path

    const uint8_t dataType = GetDataType<T>(); // dataType
    CopyToBuffer(buffer, position, &dataType);

    constexpr char no = 'n'; // isDimension
    CopyToBuffer(buffer, position, &no);

    const uint8_t dimensions = variable.m_Count.size();
    CopyToBuffer(buffer, position, &dimensions); // count

    // 27 is from 9 bytes for each: var y/n + local, var y/n + global dimension,
    // var y/n + global offset, changed for characteristic
    uint16_t dimensionsLength = 27 * dimensions;
    CopyToBuffer(buffer, position, &dimensionsLength); // length

    WriteDimensionsRecord(variable.m_Count, variable.m_Shape, variable.m_Start,
                          buffer, position);

    // CHARACTERISTICS
    WriteVariableCharacteristics(variable, stats, buffer, position);

    // Back to varLength including payload size
    // not need to remove its own size (8) from length from bpdump
    const uint64_t varLength =
        position - varLengthPosition + variable.PayLoadSize();

    size_t backPosition = varLengthPosition;
    CopyToBuffer(buffer, backPosition, &varLength);

    m_HeapBuffer.m_DataAbsolutePosition += position - varLengthPosition;
}

template <class T>
void BP1Writer::WriteVariableMetadataInIndex(
    const Variable<T> &variable,
    const Stats<typename TypeInfo<T>::ValueType> &stats, const bool isNew,
    BP1Index &index) noexcept
{
    auto &buffer = index.Buffer;

    if (isNew) // write variable header (might be shared with
               // attributes index)
    {
        buffer.insert(buffer.end(), 4, 0); // skip var length (4)
        InsertToBuffer(buffer, &stats.MemberID);
        buffer.insert(buffer.end(), 2, 0); // skip group name
        WriteNameRecord(variable.m_Name, buffer);
        buffer.insert(buffer.end(), 2, 0); // skip path

        const std::uint8_t dataType = GetDataType<T>();
        InsertToBuffer(buffer, &dataType);

        // Characteristics Sets Count in Metadata
        index.Count = 1;
        InsertToBuffer(buffer, &index.Count);
    }
    else // update characteristics sets count
    {
        if (m_Verbosity == 0)
        {
            ++index.Count;
            size_t setsCountPosition = 15 + variable.m_Name.size();
            CopyToBuffer(buffer, setsCountPosition, &index.Count);
        }
    }

    WriteVariableCharacteristics(variable, stats, buffer);
}

template <class T>
void BP1Writer::WriteBoundsRecord(const bool isScalar, const Stats<T> &stats,
                                  std::uint8_t &characteristicsCounter,
                                  std::vector<char> &buffer) noexcept
{
    if (isScalar)
    {
        // stats.min = stats.max = value, need to test
        WriteCharacteristicRecord(characteristic_value, characteristicsCounter,
                                  stats.Min, buffer);
    }
    else
    {
        if (m_Verbosity == 0) // default verbose
        {
            WriteCharacteristicRecord(
                characteristic_min, characteristicsCounter, stats.Min, buffer);

            WriteCharacteristicRecord(
                characteristic_max, characteristicsCounter, stats.Max, buffer);
        }
    }
}

template <class T>
void BP1Writer::WriteBoundsRecord(const bool singleValue, const Stats<T> &stats,
                                  std::uint8_t &characteristicsCounter,
                                  std::vector<char> &buffer,
                                  size_t &position) noexcept
{
    if (singleValue)
    {
        // stats.min = stats.max = value, need to test
        WriteCharacteristicRecord(characteristic_value, characteristicsCounter,
                                  stats.Min, buffer, position);
    }
    else
    {
        if (m_Verbosity == 0) // default min and max only
        {
            WriteCharacteristicRecord(characteristic_min,
                                      characteristicsCounter, stats.Min, buffer,
                                      position);

            WriteCharacteristicRecord(characteristic_max,
                                      characteristicsCounter, stats.Max, buffer,
                                      position);
        }
    }
}

template <class T>
void BP1Writer::WriteCharacteristicRecord(const std::uint8_t characteristicID,
                                          std::uint8_t &characteristicsCounter,
                                          const T &value,
                                          std::vector<char> &buffer) noexcept
{
    const std::uint8_t id = characteristicID;
    InsertToBuffer(buffer, &id);
    InsertToBuffer(buffer, &value);
    ++characteristicsCounter;
}

template <class T>
void BP1Writer::WriteCharacteristicRecord(const uint8_t characteristicID,
                                          uint8_t &characteristicsCounter,
                                          const T &value,
                                          std::vector<char> &buffer,
                                          size_t &position) noexcept
{
    const std::uint8_t id = characteristicID;
    CopyToBuffer(buffer, position, &id);
    CopyToBuffer(buffer, position, &value);
    ++characteristicsCounter;
}

template <class T>
void BP1Writer::WriteVariableCharacteristics(
    const Variable<T> &variable,
    const Stats<typename TypeInfo<T>::ValueType> &stats,
    std::vector<char> &buffer) noexcept
{
    // going back at the end
    const size_t characteristicsCountPosition = buffer.size();
    // skip characteristics count(1) + length (4)
    buffer.insert(buffer.end(), 5, 0);
    uint8_t characteristicsCounter = 0;

    // DIMENSIONS
    uint8_t characteristicID = characteristic_dimensions;
    InsertToBuffer(buffer, &characteristicID);
    const uint8_t dimensions = variable.m_Count.size();
    InsertToBuffer(buffer, &dimensions); // count
    const uint16_t dimensionsLength = 24 * dimensions;
    InsertToBuffer(buffer, &dimensionsLength); // length
    WriteDimensionsRecord(variable.m_Count, variable.m_Shape, variable.m_Start,
                          buffer);
    ++characteristicsCounter;

    WriteBoundsRecord(variable.m_SingleValue, stats, characteristicsCounter,
                      buffer);

    WriteCharacteristicRecord(characteristic_time_index, characteristicsCounter,
                              stats.TimeIndex, buffer);

    WriteCharacteristicRecord(characteristic_offset, characteristicsCounter,
                              stats.Offset, buffer);

    WriteCharacteristicRecord(characteristic_payload_offset,
                              characteristicsCounter, stats.PayloadOffset,
                              buffer);
    // END OF CHARACTERISTICS

    // Back to characteristics count and length
    size_t backPosition = characteristicsCountPosition;
    CopyToBuffer(buffer, backPosition, &characteristicsCounter); // count (1)

    // remove its own length (4) + characteristic counter (1)
    const uint32_t characteristicsLength =
        buffer.size() - characteristicsCountPosition - 4 - 1;

    CopyToBuffer(buffer, backPosition, &characteristicsLength); // length
}

template <class T>
void BP1Writer::WriteVariableCharacteristics(
    const Variable<T> &variable,
    const Stats<typename TypeInfo<T>::ValueType> &stats,
    std::vector<char> &buffer, size_t &position) noexcept
{
    // going back at the end
    const size_t characteristicsCountPosition = position;
    // skip characteristics count(1) + length (4)
    position += 5;
    uint8_t characteristicsCounter = 0;

    // DIMENSIONS
    uint8_t characteristicID = characteristic_dimensions;
    CopyToBuffer(buffer, position, &characteristicID);

    const uint8_t dimensions = variable.m_Count.size();
    CopyToBuffer(buffer, position, &dimensions); // count
    const uint16_t dimensionsLength = 24 * dimensions;
    CopyToBuffer(buffer, position, &dimensionsLength); // length
    WriteDimensionsRecord(variable.m_Count, variable.m_Shape, variable.m_Start,
                          buffer, position, true); // isCharacteristic = true
    ++characteristicsCounter;

    // VALUE for SCALAR or STAT min, max for ARRAY
    WriteBoundsRecord(variable.m_SingleValue, stats, characteristicsCounter,
                      buffer, position);
    // END OF CHARACTERISTICS

    // Back to characteristics count and length
    size_t backPosition = characteristicsCountPosition;
    CopyToBuffer(buffer, backPosition, &characteristicsCounter);

    // remove its own length (4) + characteristic counter (1)
    const uint32_t characteristicsLength =
        position - characteristicsCountPosition - 4 - 1;
    CopyToBuffer(buffer, backPosition, &characteristicsLength);
}

template <class T>
void BP1Writer::WriteAttributeInData(const Attribute<T> &attribute,
                                     Stats<T> &stats) noexcept
{
    auto &buffer = m_HeapBuffer.m_Data;
    auto &position = m_HeapBuffer.m_DataPosition;

    // for writing length at the end
    const size_t attributeLengthPosition = position;
    position += 4; // skip attribute length (4)

    // CopyToBuffer(buffer, position, &memberID);

    WriteNameRecord(attribute.m_Name, buffer, position);
    position += 2; // skip path

    constexpr char isVariableRelated = 'n';
    CopyToBuffer(buffer, position, &isVariableRelated);

    const int8_t dataType = GetDataType<T>();
    // single string, only possible using vector
    if (dataType == type_string)
    {
        const uint32_t dataSize =
            static_cast<uint32_t>(attribute.m_Data->size());

        if (dataSize == 1) // string string
        {
            CopyToBuffer(buffer, position, &dataType); // type_string

            const std::string singleString(
                reinterpret_cast<const std::string>(attribute.m_Data->front()));

            const uint32_t length = singleString.length();
            CopyToBuffer(buffer, position, &length);
            CopyToBuffer(buffer, position, singleString.c_str(), length);
        }
        else if (dataSize > 1) // string array
        {
            const int8_t stringArrayDataType =
                static_cast<int8_t>(type_string_array);
            CopyToBuffer(buffer, position,
                         &stringArrayDataType); // type_string_array

            CopyToBuffer(buffer, position, &dataSize);

            std::vector<std::string> stringsData =
                *reinterpret_cast<std::vector<std::string> *>(attribute.m_Data);

            for (const auto &data : stringsData)
            {
                const std::string singleString(
                    reinterpret_cast<std::string>(data) + '\0');

                const uint32_t length = singleString.length();
                CopyToBuffer(buffer, position, &length);
                CopyToBuffer(buffer, position, singleString.c_str(), length);
            }
        }
    }
    else // non-string data array
    {
        if (attribute.m_IsArray) // non-string
        {
            const uint32_t length =
                static_cast<uint32_t>(attribute.m_Elements * sizeof(T));
            CopyToBuffer(buffer, position, &length);
            CopyToBuffer(buffer, position, attribute.m_ArrayData, length);
        }
        else
        {
            const uint32_t length =
                static_cast<uint32_t>(attribute.m_Data->size() * sizeof(T));
            CopyToBuffer(buffer, position, &length);
            CopyToBuffer(buffer, position, attribute.m_Data->data(), length);
        }
    }

    size_t backPosition = attributeLengthPosition;
    // back to write length TODO check with bpdump to keep or remove 4
    const uint32_t attributeLength = position - attributeLengthPosition - 4;
    CopyToBuffer(buffer, backPosition, &attributeLength); // length
    // update absolute position
    m_HeapBuffer.m_DataAbsolutePosition = position - attributeLengthPosition;
}

template <class T>
void BP1Writer::WriteAttributeInIndex(const Attribute<T> &attribute,
                                      const Stats<T> &stats,
                                      BP1Index &index) noexcept
{
    auto &buffer = index.Buffer;

    buffer.insert(buffer.end(), 4, 0); // skip var length (4)
    InsertToBuffer(buffer, &stats.MemberID);
    buffer.insert(buffer.end(), 2, 0); // skip group name
    WriteNameRecord(attribute.m_Name, buffer);
    buffer.insert(buffer.end(), 2, 0); // skip path

    const std::uint8_t dataType = GetDataType<T>();
    InsertToBuffer(buffer, &dataType);

    index.Count = 1;
    InsertToBuffer(buffer, &index.Count);

    // TODO finish characteristics in index
}

} // end namespace format
} // end namespace adios

#endif // ADIOS2_UTILITIES_FORMAT_BP1_BP1WRITER_TCC_
