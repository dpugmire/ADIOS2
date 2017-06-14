/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * BP1Base.h  base class for BP1Writer and BP1Reader
 *
 *  Created on: Feb 2, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_FORMAT_BP1_BP1BASE_H_
#define ADIOS2_TOOLKIT_FORMAT_BP1_BP1BASE_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <memory> //std::shared_ptr
#include <vector>
/// \endcond

#include "adios2/ADIOSConfig.h"
#include "adios2/ADIOSMPICommOnly.h"
#include "adios2/ADIOSMacros.h"
#include "adios2/ADIOSTypes.h"
#include "adios2/core/Variable.h"
#include "adios2/toolkit/capsule/heap/STLVector.h"
#include "adios2/toolkit/format/bp1/BP1Aggregator.h"
#include "adios2/toolkit/format/bp1/BP1Structs.h"
#include "adios2/toolkit/profiling/iochrono/IOChrono.h"

namespace adios
{
namespace format
{

/**
 * Base class for BP1Writer and BP1Reader format
 */
class BP1Base
{

public:
    /** statistics verbosity, only 0 is supported */
    unsigned int m_Verbosity = 0;

    /** contains data buffer and position */
    capsule::STLVector m_HeapBuffer;

    /** memory growth factor,s set by the user */
    float m_GrowthFactor = DefaultBufferGrowthFactor;

    /** max buffer size, set by the user */
    size_t m_MaxBufferSize = DefaultMaxBufferSize;

    /** contains bp1 format metadata indices*/
    BP1MetadataSet m_MetadataSet;

    /** object that takes care of all MPI aggregation tasks */
    BP1Aggregator m_BP1Aggregator;

    /** true: Close was called, Engine will call this many times for different
     * transports */
    bool m_IsClosed = false;

    /** buffering and MPI aggregation profiling info, set by user */
    profiling::IOChrono m_Profiler;

    /**
     * Unique constructor
     * @param mpiComm for m_BP1Aggregator
     * @param debugMode true: exceptions checks
     */
    BP1Base(MPI_Comm mpiComm, const bool debugMode);

    virtual ~BP1Base() = default;

    void InitParameters(const Params &parameters);

    /**
     * Vector version of BPBaseName
     * @param names
     * @return vector of base (name.bp) names
     */
    std::vector<std::string>
    GetBPBaseNames(const std::vector<std::string> &names) const noexcept;

    /**
     * Vector version of GetBPName
     * @param names
     * @return
     */
    std::vector<std::string>
    GetBPNames(const std::vector<std::string> &baseNames) const noexcept;

    /** Return type of the CheckAllocation function. */
    enum class ResizeResult
    {
        Failure,   //!< FAILURE, caught a std::bad_alloc
        Unchanged, //!< UNCHANGED, no need to resize (sufficient capacity)
        Success,   //!< SUCCESS, resize was successful
        Flush      //!< FLUSH, need to flush to transports for current variable
    };

    /**
     * @param variable
     * @return
     * -1: allocation failed,
     *  0: no allocation needed,
     *  1: reallocation is sucessful
     *  2: need a transport flush
     */
    template <class T>
    ResizeResult ResizeBuffer(const Variable<T> &variable);

protected:
    /** might be used in large payload copies to buffer */
    unsigned int m_Threads = 1;
    const bool m_DebugMode = false;

    /** method type for file I/O */
    enum IO_METHOD
    {
        METHOD_UNKNOWN = -2,
        METHOD_NULL = -1,
        METHOD_MPI = 0,
        METHOD_DATATAP = 1 // OBSOLETE
        ,
        METHOD_POSIX = 2,
        METHOD_DATASPACES = 3,
        METHOD_VTK = 4 // non-existent
        ,
        METHOD_POSIX_ASCII = 5 // non-existent
        ,
        METHOD_MPI_CIO = 6 // OBSOLETE
        ,
        METHOD_PHDF5 = 7,
        METHOD_PROVENANCE = 8 // OBSOLETE
        ,
        METHOD_MPI_STRIPE = 9 // OBSOLETE
        ,
        METHOD_MPI_LUSTRE = 10,
        METHOD_MPI_STAGGER = 11 // OBSOLETE
        ,
        METHOD_MPI_AGG = 12 // OBSOLETE
        ,
        METHOD_ADAPTIVE = 13 // OBSOLETE
        ,
        METHOD_POSIX1 = 14 // OBSOLETE
        ,
        METHOD_NC4 = 15,
        METHOD_MPI_AMR = 16,
        METHOD_MPI_AMR1 = 17 // OBSOLETE
        ,
        METHOD_FLEXPATH = 18,
        METHOD_NSSI_STAGING = 19,
        METHOD_NSSI_FILTER = 20,
        METHOD_DIMES = 21,
        METHOD_VAR_MERGE = 22,
        METHOD_MPI_BGQ = 23,
        METHOD_ICEE = 24,
        METHOD_COUNT = 25,
        METHOD_FSTREAM = 26,
        METHOD_FILE = 27,
        METHOD_ZMQ = 28,
        METHOD_MDTM = 29
    };

    /**
     * DataTypes mapping in BP Format
     */
    enum DataTypes
    {
        type_unknown = -1, //!< type_unknown
        type_byte = 0,     //!< type_byte
        type_short = 1,    //!< type_short
        type_integer = 2,  //!< type_integer
        type_long = 4,     //!< type_long

        type_unsigned_byte = 50,    //!< type_unsigned_byte
        type_unsigned_short = 51,   //!< type_unsigned_short
        type_unsigned_integer = 52, //!< type_unsigned_integer
        type_unsigned_long = 54,    //!< type_unsigned_long

        type_real = 5,        //!< type_real or float
        type_double = 6,      //!< type_double
        type_long_double = 7, //!< type_long_double

        type_string = 9,              //!< type_string
        type_complex = 10,            //!< type_complex
        type_double_complex = 11,     //!< type_double_complex
        type_string_array = 12,       //!< type_string_array
        type_long_double_complex = 13 //!< type_double_complex
    };

    /**
     * Characteristic ID in variable metadata
     */
    enum VariableCharacteristicID
    {
        characteristic_value = 0,      //!< characteristic_value
        characteristic_min = 1,        //!< Used to read in older bp file format
        characteristic_max = 2,        //!< Used to read in older bp file format
        characteristic_offset = 3,     //!< characteristic_offset
        characteristic_dimensions = 4, //!< characteristic_dimensions
        characteristic_var_id = 5,     //!< characteristic_var_id
        characteristic_payload_offset = 6, //!< characteristic_payload_offset
        characteristic_file_index = 7,     //!< characteristic_file_index
        characteristic_time_index = 8,     //!< characteristic_time_index
        characteristic_bitmap = 9,         //!< characteristic_bitmap
        characteristic_stat = 10,          //!< characteristic_stat
        characteristic_transform_type = 11 //!< characteristic_transform_type
    };

    /** Define statistics type for characteristic ID = 10 in bp1 format */
    enum VariableStatistics
    {
        statistic_min = 0,
        statistic_max = 1,
        statistic_cnt = 2,
        statistic_sum = 3,
        statistic_sum_square = 4,
        statistic_hist = 5,
        statistic_finite = 6
    };

    template <class T>
    struct Stats
    {
        T Min;
        T Max;
        uint64_t Offset;
        uint64_t PayloadOffset;
        uint32_t TimeIndex;
        uint32_t MemberID;

        //      unsigned long int count;
        //      long double sum;
        //      long double sumSquare;
        // unsigned long int histogram
        // bool finite??
    };

    /** profile=on (default) generate profiling.log
         *  profile=off */
    void InitParameterProfile(const std::string value);

    /** profile_units=s (default) (mus, ms, s,m,h) from ADIOSTypes.h TimeUnit */
    void InitParameterProfileUnits(const std::string value);

    /** growth_factor=1.5 (default), must be > 1.0 */
    void InitParameterBufferGrowth(const std::string value);

    /** set initial buffer size */
    void InitParameterInitBufferSize(const std::string value);

    /** unlimited (default), set max buffer size in Gb or Mb
     *  max_buffer_size=100Mb or  max_buffer_size=1Gb */
    void InitParameterMaxBufferSize(const std::string value);

    /** verbose file level=0 (default) */
    void InitParameterVerbose(const std::string value);

    /**
     * Returns data type index from enum Datatypes
     * @param variable input variable
     * @return data type
     */
    template <class T>
    int8_t GetDataType() const noexcept;

    std::vector<uint8_t>
    GetTransportIDs(const std::vector<std::string> &transportsTypes) const
        noexcept;

    /**
     * Calculates the Process Index size in bytes according to the BP
     * format,
     * including list of method with no parameters (for now)
     * @param name
     * @param timeStepName
     * @param transportsSize
     * @return size of pg index
     */
    size_t GetProcessGroupIndexSize(const std::string name,
                                    const std::string timeStepName,
                                    const size_t transportsSize) const noexcept;

    /**
     * Returns the estimated variable index size
     * @param variable
     */
    template <class T>
    size_t GetVariableIndexSize(const Variable<T> &variable) const noexcept;
};

#define declare_template_instantiation(T)                                      \
    extern template BP1Base::ResizeResult BP1Base::ResizeBuffer(               \
        const Variable<T> &variable);

ADIOS2_FOREACH_TYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation

} // end namespace format
} // end namespace adios

#endif /* ADIOS2_TOOLKIT_FORMAT_BP1_BP1BASE_H_ */
