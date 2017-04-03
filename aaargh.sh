#!/bin/bash

BASEDIR=$(readlink -f $(dirname ${BASH_SOURCE}))
cd ${BASEDIR}

echo "Updating scripts..."
git pull --ff-only
cd ${HOME}
mkdir -p ${BASEDIR}/../Logs


for COMP_ID in Intel GCC
do
  for MPI_ID in OpenMPI MPICH MVAPICH2 NoMPI
  do
    module purge
    module load ohpc

    case ${MPI_ID} in
      OpenMPI)  ;;
      MPICH)    module swap openmpi mpich ;;
      MVAPICH2) module swap openmpi mvapich2 ;;
      *) ;;
    esac

    case ${COMP_ID} in
      Intel)
        COMP_VER=$(module list | sed -n 's|.* intel\/\([^ ]*\) .*|\1|p')
        export CC=$(which icc) CXX=$(which icpc) FC=$(which ifort)
        ;;
      GCC)
        module -q swap intel gnu
        COMP_VER=$(module list | sed -n 's|.* gnu\/\([^ ]*\) .*|\1|p')
        export CC=$(which gcc) CXX=$(which g++) FC=$(which gfortran)
        ;;
      *) ;;
    esac

    echo "Building ${COMP_ID}-${COMP_VER} ${MPI_ID}"

    LOGBASE=${BASEDIR}/../Logs/${COMP_ID}-${COMP_VER}_${MPI_ID}
    ctest -DCOMP_ID=${COMP_ID} -DCOMP_VER=${COMP_VER} -DMPI_ID=${MPI_ID} \
      -S ${BASEDIR}/aaargh.cmake -VV 1>${LOGBASE}.out 2>${LOGBASE}.err
  done
done

module purge
COMP_ID=GCC
MPI_ID=NoMPI
export CC=gcc CXX=g++ FC=gfortran

COMP_VER=$(gcc --version | head -1 | awk '{print $3}')
echo "Building ${COMP_ID}-${COMP_VER} ${MPI_ID}"
LOGBASE=${BASEDIR}/../Logs/${COMP_ID}-${COMP_VER}_${MPI_ID}
ctest \
  -DCOMP_ID=${COMP_ID} -DCOMP_VER=${COMP_VER} -DMPI_ID=${MPI_ID} \
  -S ${BASEDIR}/aaargh.cmake -VV 1>${LOGBASE}.out 2>${LOGBASE}.err

COMP_VER=$(scl enable devtoolset-3 -- gcc --version | head -1 | awk '{print $3}')
echo "Building ${COMP_ID}-${COMP_VER} ${MPI_ID}"
LOGBASE=${BASEDIR}/../Logs/${COMP_ID}-${COMP_VER}_${MPI_ID}
scl enable devtoolset-3 -- ctest \
  -DCOMP_ID=${COMP_ID} -DCOMP_VER=${COMP_VER} -DMPI_ID=${MPI_ID} \
  -S ${BASEDIR}/aaargh.cmake -VV 1>${LOGBASE}.out 2>${LOGBASE}.err

COMP_VER=$(scl enable devtoolset-4 -- gcc --version | head -1 | awk '{print $3}')
echo "Building ${COMP_ID}-${COMP_VER} ${MPI_ID}"
LOGBASE=${BASEDIR}/../Logs/${COMP_ID}-${COMP_VER}_${MPI_ID}
scl enable devtoolset-4 -- ctest \
  -DCOMP_ID=${COMP_ID} -DCOMP_VER=${COMP_VER} -DMPI_ID=${MPI_ID} \
  -S ${BASEDIR}/aaargh.cmake -VV 1>${LOGBASE}.out 2>${LOGBASE}.err

COMP_VER=$(scl enable devtoolset-6 -- gcc --version | head -1 | awk '{print $3}')
echo "Building ${COMP_ID}-${COMP_VER} ${MPI_ID}"
LOGBASE=${BASEDIR}/../Logs/${COMP_ID}-${COMP_VER}_${MPI_ID}
scl enable devtoolset-6 -- ctest \
  -DCOMP_ID=${COMP_ID} -DCOMP_VER=${COMP_VER} -DMPI_ID=${MPI_ID} \
  -S ${BASEDIR}/aaargh.cmake -VV 1>${LOGBASE}.out 2>${LOGBASE}.err

echo "Building ${COMP_ID}-${COMP_VER} ${MPI_ID} ClangTidy"
LOGBASE=${BASEDIR}/../Logs/${COMP_ID}-${COMP_VER}_${MPI_ID}_ClangTidy
ctest -DCOMP_ID=${COMP_ID} -DCOMP_VER=${COMP_VER} -DMPI_ID=${MPI_ID} \
  -S ${BASEDIR}/aaargh_clang_tidy.cmake -VV 1>${LOGBASE}.out 2>${LOGBASE}.err

export CFLAGS="-fprofile-arcs -ftest-coverage"
export CXXFLAGS="-fprofile-arcs -ftest-coverage"
COMP_VER=$(gcc --version | head -1 | awk '{print $3}')
echo "Building ${COMP_ID}-${COMP_VER} ${MPI_ID} Debug"
LOGBASE=${BASEDIR}/../Logs/${COMP_ID}-${COMP_VER}_${MPI_ID}
ctest \
  -DCOMP_ID=${COMP_ID} -DCOMP_VER=${COMP_VER} -DMPI_ID=${MPI_ID} \
  -S ${BASEDIR}/aaargh_memcheck_coverage.cmake -VV 1>${LOGBASE}.out 2>${LOGBASE}.err
