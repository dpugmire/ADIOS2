#------------------------------------------------------------------------------#
# Distributed under the OSI-approved Apache License, Version 2.0.  See
# accompanying file Copyright.txt for details.
#------------------------------------------------------------------------------#
  
add_library(adios2
  core/Attribute.cpp core/Attribute.tcc
  core/AttributeBase.cpp
  core/ADIOS.cpp
  core/Engine.cpp
  core/IO.cpp core/IO.tcc
  core/Selection.cpp
  core/SelectionBoundingBox.cpp
  core/SelectionPoints.cpp
  core/Transform.cpp core/Transform.tcc
  core/Variable.cpp core/Variable.tcc
  core/VariableBase.cpp
  core/VariableCompound.cpp core/VariableCompound.tcc
  
#helper
  helper/adiosMath.cpp
  helper/adiosMPIFunctions.cpp
  helper/adiosString.cpp
  helper/adiosSystem.cpp
  helper/adiosType.cpp
  helper/adiosXML.cpp
  
#  engine/bp/BPFileReader.cpp
  engine/bp/BPFileWriter.cpp engine/bp/BPFileWriter.tcc

  toolkit/capsule/Capsule.cpp
  toolkit/capsule/heap/STLVector.cpp 
 
  toolkit/format/bp1/BP1Base.cpp toolkit/format/bp1/BP1Base.tcc
  toolkit/format/bp1/BP1Aggregator.cpp
  toolkit/format/bp1/BP1Writer.cpp toolkit/format/bp1/BP1Writer.tcc

  toolkit/profiling/iochrono/Timer.cpp

  toolkit/transport/Transport.cpp
  toolkit/transport/file/FilePointer.cpp
  toolkit/transport/file/FileStream.cpp
  
  toolkit/transportman/TransportMan.cpp  
)
target_include_directories(adios2
  PUBLIC
    $<BUILD_INTERFACE:${ADIOS2_SOURCE_DIR}/source>
    $<BUILD_INTERFACE:${ADIOS2_BINARY_DIR}/source>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(adios2 PRIVATE adios2sys pugixml)

find_package(Threads REQUIRED)
target_link_libraries(adios2 PUBLIC ${CMAKE_THREAD_LIBS_INIT})

if(UNIX)
  target_sources(adios2 PRIVATE toolkit/transport/file/FileDescriptor.cpp)
endif()

if(ADIOS2_HAVE_SysVShMem)
  target_sources(adios2 PRIVATE toolkit/capsule/shmem/ShmSystemV.cpp)
endif()

if(ADIOS2_HAVE_DataMan)
  target_sources(adios2 PRIVATE
    engine/dataman/DataManReader.cpp
    engine/dataman/DataManWriter.cpp
  )
  target_link_libraries(adios2 PRIVATE dataman)
endif()
  
if(ADIOS2_HAVE_BZip2)
  find_package(BZip2 REQUIRED)
  target_sources(adios2 PRIVATE transform/compress/CompressBZip2.cpp)
  target_link_libraries(adios2 PRIVATE BZip2::BZip2)
endif()

if(ADIOS2_HAVE_ZFP)
  find_package(ZFP REQUIRED)
  target_sources(adios2 PRIVATE transform/compress/CompressZfp.cpp)
  target_link_libraries(adios2 PRIVATE zfp::zfp)
endif()

if(ADIOS2_HAVE_VTKm)
  find_package (VTKm REQUIRED OPTIONAL_COMPONENTS Serial OpenGL Rendering TBB CUDA)
  target_sources(adios2 PRIVATE toolkit/analytics/vis/VisVTKm.cpp)
  target_link_libraries(adios2 PRIVATE ${VTKm_LIBRARIES})
endif()

if(ADIOS2_HAVE_MPI)
  find_package(MPI COMPONENTS C REQUIRED)
  target_include_directories(adios2 PUBLIC ${MPI_C_INCLUDE_PATH})
  target_link_libraries(adios2 PUBLIC ${MPI_C_LIBRARIES})
else()
  target_sources(adios2 PRIVATE mpidummy.cpp)
endif()

if(ADIOS2_HAVE_ADIOS1)
  if(ADIOS2_HAVE_MPI)
    find_package(ADIOS1 1.12.0 REQUIRED)
  else()
    find_package(ADIOS1 1.12.0 COMPONENTS sequential REQUIRED)
  endif()

  target_sources(adios2 PRIVATE
    engine/adios1/ADIOS1Reader.cpp
    engine/adios1/ADIOS1Writer.cpp
    
    toolkit/interop/adios1/ADIOS1Common.cpp toolkit/interop/adios1/ADIOS1Common.tcc
    
  )
  target_link_libraries(adios2 PRIVATE adios1::adios)
endif()

if(ADIOS2_HAVE_HDF5)
  find_package(HDF5 REQUIRED)
  if(ADIOS2_HAVE_MPI AND (NOT HDF5_IS_PARALLEL))
    message(FATAL_ERROR
      "A sequential version of HDF5 was detected but the parallel version "
      "of ADIOS is being built, which requires a parallel HDF5."
    )
  elseif((NOT ADIOS2_HAVE_MPI) AND HDF5_IS_PARALLEL)
    message(FATAL_ERROR
      "A parallel version of HDF5 was detected but the sequential version "
      "of ADIOS is being built, which requires a sequential HDF5."
    )
  endif()
  if(HDF5_C_INCLUDE_DIRS)
    target_include_directories(adios2 PRIVATE ${HDF5_C_INCLUDE_DIRS})
  else()
    target_include_directories(adios2 PRIVATE ${HDF5_INCLUDE_DIRS})
  endif()
  target_sources(adios2 PRIVATE
    engine/hdf5/HDF5ReaderP.cpp
    engine/hdf5/HDF5WriterP.cpp
    toolkit/interop/hdf5/HDF5Common.cpp toolkit/interop/hdf5/HDF5Common.tcc
  )
  target_link_libraries(adios2 PRIVATE ${HDF5_C_LIBRARIES})
endif()

# Set library version information
set_target_properties(adios2 PROPERTIES
  VERSION ${ADIOS2_VERSION}
  SOVERSION ${ADIOS2_VERSION_MAJOR}
)

install(FILES ADIOSMacros.h ADIOSTypes.h ADIOSMPICommOnly.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/adios2
)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/core
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/adios2
  FILES_MATCHING REGEX "[^/]*\.(h|inl)$"
)

# Library installation
install(TARGETS adios2 EXPORT adios2Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
