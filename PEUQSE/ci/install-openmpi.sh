# Variables that should be set:
# MPI_BUILD_DIR=/path/to/install/openmpi
# MPI_VERSION="3.1" or some other supported version (see open-mpi.org)
# MPI_FULL_VERSION="3.1.3" or some other supported version.
# C_COMPILER="path to C compiler to use for mpicc"
# CXX_COMPILER="path to C++ compiler to use for mpic++"
#
# After sourcing this script, the following variables will be set/altered:
# MPI_CC : Path to MPI C compiler
# MPI_CXX : Path to MPI C++ compiler
# MPI_EXEC : Path to the corresponding mpiexec launcher

CC_OLD=$CC
CXX_OLD=$CXX
export CC=$C_COMPILER
export CXX=$CXX_COMPILER

# Check to see if we have a cached build from earlier.
if [ -f "$MPI_BUILD_DIR/bin/mpiexec" ] && [ -f "openmpi-$MPI_FULL_VERSION/config.log" ]; then
    echo "Using cached OpenMPI"
else
    echo "Downloading OpenMPI Source"
    wget https://download.open-mpi.org/release/open-mpi/v$MPI_VERSION/openmpi-$MPI_FULL_VERSION.tar.gz
    tar zxf openmpi-$MPI_FULL_VERSION.tar.gz
    echo "Configuring and building OpenMPI"
    mkdir -p $MPI_BUILD_DIR
    cd openmpi-$MPI_FULL_VERSION
    export CC=$C_COMPILER
    export CXX=$CXX_COMPILER
    ./configure --prefix=$MPI_BUILD_DIR &> config.log
    make -j3  # Produce output so that Travis doesn't halt after 10 min.
    make install &> make-install.log
    cd ..
fi
export MPI_CC=$MPI_BUILD_DIR/bin/mpicc
export MPI_CXX=$MPI_BUILD_DIR/bin/mpic++
export MPI_EXEC=$MPI_BUILD_DIR/bin/mpiexec
export CC=$CC_OLD
export CXX=$CXX_OLD
