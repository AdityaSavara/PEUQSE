#!/bin/bash
#Following: https://stackoverflow.com/questions/50687222/parallel-group-setup-mpi4py-openmdao-2-2-x
wget https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.0.tar.gz 
tar -xzf openmpi-3.1.0.tar.gz 
cd openmpi-*
./configure --prefix="/home/$USER/.openmpi"
make
sudo make install
echo export PATH="$PATH:/home/$USER/.openmpi/bin" >> /home/$USER/.bashrc
echo export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/$USER/.openmpi/lib/" >> /home/$USER/.bashrc