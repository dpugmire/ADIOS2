/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * TransportMan.h : manages a vector of transports
 *
 *  Created on: May 23, 2017
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_TOOLKIT_TRANSPORT_TRANSPORTMANAGER_H_
#define ADIOS2_TOOLKIT_TRANSPORT_TRANSPORTMANAGER_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <memory> //std::shared_ptr
#include <string>
#include <unordered_map>
#include <vector>
/// \endcond

#include "adios2/ADIOSMPICommOnly.h"
#include "adios2/toolkit/transport/Transport.h"

namespace adios
{
namespace transportman
{

class TransportMan
{

public:
    bool m_GetCollectiveMetadata = false;

    /**
     * Unique base constructor
     * @param mpiComm
     * @param debugMode
     */
    TransportMan(MPI_Comm mpiComm, const bool debugMode);

    virtual ~TransportMan() = default;

    void SetParameters(const std::vector<Params> &parametersVector);

    /**
     * Gets each transport base name from either baseName at Open or name key in
     * parameters
     * Checks if transport name rules IO AddTransport have unique names for
     * every type (for now)
     * @param baseName from Open
     * @param parameters from IO TransportsParameters (from AddTransport
     * function)
     * @return transport base names
     */
    std::vector<std::string>
    GetFilesBaseNames(const std::string &baseName,
                      const std::vector<Params> &parametersVector) const;

    /**
     * Open transport files from IO AddTransport
     * @param baseNames passed from Open( name )
     * @param names actual filenames (from BP)
     * @param openMode
     * @param parametersVector from IO
     * @param profile
     */
    void OpenFiles(const std::vector<std::string> &baseNames,
                   const std::vector<std::string> &names,
                   const OpenMode openMode,
                   const std::vector<Params> &parametersVector,
                   const bool profile);

    /**
     * Checks if index is in range
     * @param index input to be checked against m_Transports range or it's -1
     * @param true: in range, false: out of range
     */
    void CheckTransportIndex(const int index) const;

    /**
     * m_Type from m_Transports based on derived classes of Transport
     * @return m_Type for each transport in m_Transports (e.g.
     * {FileDescriptor,
     * FilePointer} )
     */
    std::vector<std::string> GetTransportsTypes() noexcept;

    /** Returns a vector of pointer references (not owning the memory) to
     * m_Transports.m_Profiler */
    std::vector<profiling::IOChrono *> GetTransportsProfilers() noexcept;

    /**
     * Write to file transports
     * @param transportIndex
     * @param buffer
     * @param size
     */
    void WriteFiles(const char *buffer, const size_t size,
                    const int transportIndex = -1);

    /**
     * Close file or files depending on transport index. Throws an exception
     * if transport is not a file when transportIndex > -1.
     * @param transportIndex -1: all transports, otherwise index in m_Transports
     */
    void CloseFiles(const int transportIndex = -1);

    /** Checks if all transports are closed */
    bool AllTransportsClosed() const noexcept;

protected:
    MPI_Comm m_MPIComm;
    const bool m_DebugMode = false;

    /** contains all transports from IO AddTransport
         * <pre>
         * key : unique id from IO AddTransport
         * value : obejct derived from Transport.h class
         * </pre>
         */
    std::vector<std::shared_ptr<Transport>> m_Transports;

    /**
     * Replaces what vector<bool> would have been
     */
    enum class CollectiveMetadata
    {
        On, //!< On
        Off //!< Off (default)
    };

    /** Keep track if transport requires collective metadata */
    std::vector<CollectiveMetadata> m_CollectiveMetadata;

    /** Called from SetParameters */
    void InitCollectiveMetadata(const Params &parameters);

    /** Called from OpenFileTransports */
    void OpenFileTransport(const std::string &fileName, const OpenMode openMode,
                           const Params &parameters, const bool profile);
};

} // end namespace transport
} // end namespace adios

#endif /* ADIOS2_TOOLKIT_TRANSPORT_TRANSPORTMANAGER_H_ */
