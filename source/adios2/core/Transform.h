/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * Transform.h : Base class for all transforms under adios2/transform
 *
 *  Created on: Oct 17, 2016
 *      Author: William F Godoy godoywf@ornl.gov
 */

#ifndef ADIOS2_CORE_TRANSFORM_H_
#define ADIOS2_CORE_TRANSFORM_H_

/// \cond EXCLUDE_FROM_DOXYGEN
#include <string>
#include <vector>
/// \endcond

#include "adios2/ADIOSTypes.h"
#include <iostream>

namespace adios2
{

/** Base class that defines data variable transformations implemented under
 * adios2/transform */
class Transform
{

public:
    /** From derived class */
    const std::string m_Library;
    const Params m_Parameters;

    /**
     * Unique base class constructor
     * @param method bzip2, zfp
     * @param debugMode true: extra exceptions checks
     */
    Transform(const std::string library, const bool debugMode);

    Transform(const std::string library, const Params &parameters,
              const bool debugMode);

    virtual ~Transform() = default;

    /**
     * Returns a conservative buffer size to hold input data for classes
     * @param sizeIn size of input data to be compressed in bytes
     * @return recommended allocation for output buffer
     */
    virtual size_t BufferMaxSize(const size_t sizeIn) const;

    /**
     * Used by Zfp
     * Returns a conservative buffer size to hold input data for classes
     * @param dataIn
     * @param dimensions
     * @return recommended allocation for output buffer in bytes
     */
    template <class T>
    size_t BufferMaxSize(const T *dataIn, const Dims &dimensions,
                         const Params &params) const;

    /**
     * BZip2 and Zfp common call
     * @param dataIn
     * @param dimensions
     * @param elementSize
     * @param type
     * @param bufferOut
     * @param parameters
     * @return size of compressed buffer
     */
    virtual size_t Compress(const void *dataIn, const Dims &dimensions,
                            const size_t elementSize, const std::string type,
                            void *bufferOut,
                            const Params &parameters = Params()) const;

    virtual size_t Decompress(const void *bufferIn, const size_t sizeIn,
                              void *dataOut, const size_t sizeOut) const;

    /**
     * Zfp signature
     * @param bufferIn
     * @param sizeIn
     * @param dataOut
     * @param dimensions
     * @param type
     * @return
     */
    virtual size_t Decompress(const void *bufferIn, const size_t sizeIn,
                              void *dataOut, const Dims &dimensions,
                              const std::string type,
                              const Params &parameters) const;

    virtual bool Render1DStructured(const void *field,
                                    const size_t fieldElements, const double x0,
                                    double deltaX)
    {
        std::cout << __FILE__ << " " << __LINE__ << std::endl;
        return false;
    }

protected:
    /** true: extra exception checks, false: skip exception checks */
    const bool m_DebugMode = false;

    /**
     * Used by CompressZfp
     * Returns a conservative buffer size to hold input data for classes
     * @param dataIn
     * @param dimensions
     * @param type
     * @return conservative buffer size for allocation
     */
    virtual size_t DoBufferMaxSize(const void *dataIn, const Dims &dimensions,
                                   const std::string type,
                                   const Params &parameters) const;
};

} // end namespace adios2

#endif /* ADIOS2_CORE_TRANSFORM_H_ */
