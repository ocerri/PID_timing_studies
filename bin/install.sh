#!/bin/bash
start_dir=$PWD
if [ ! -d "$1" ]; then
  mkdir $1
fi
cd $1

#Getting root 6.10.08
if [ -d "/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt" ]; then
  . /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt/setup.sh
  . /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.10.08/x86_64-centos7-gcc48-opt/root/bin/thisroot.sh
else
  read -p "Have you set the desired root? [y/n]" aux_rpl
  if [ aux_rpl=="y" ]; then
    echo "Root location:"
    echo $ROOTSYS
  else
    echo "Installing ROOT 6.10.08"
    if [ ! -d "root-6.10.08" ]; then
      wget https://root.cern.ch/download/root_v6.10.08.source.tar.gz
      tar -zxf root_v6.10.08.source.tar.gz
      rm -rf root_v6.10.08.source.tar.gz
      mv root-6.10.08 source_root-6.10.08
      mkdir root-6.10.08
      cd root-6.10.08
      cmake ../source_root-6.10.08/ -Dall=ON
      cmake --build . -- -j8
      source ./bin/thisroot.sh
      cd ..
      rm -rf source_root-6.10.08
    fi
  fi
fi

echo "Installing HepMC"
mkdir hepmc_20609
cd hepmc_20609
wget http://hepmc.web.cern.ch/hepmc/releases/hepmc2.06.09.tgz
tar zxf hepmc2.06.09.tgz
rm -rf hepmc2.06.09.tgz
mv hepmc2.06.09 source_hepmc20609
mkdir build
cd build
cmake ../source_hepmc20609 -DCMAKE_INSTALL_PREFIX=./ -Dmomentum:STRING=GEV -Dlength:STRING=MM
make -j8
make install
cd ..
echo "#"'!'"/bin/sh" > user_setup.sh
echo "export HEPMC2=$PWD/build" >> user_setup.sh
echo "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\$HEPMC2/lib" >> user_setup.sh
source ./user_setup.sh
cd ..

echo "Installing Pythia 8.2.3"
if [ ! -d "pythia8230" ]; then
  wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8230.tgz
  tar xvfz pythia8230.tgz
  rm -rf pythia8230.tgz
fi
cd pythia8230
./configure --with-hepmc2=$HEPMC2
gmake -j8
echo "#"'!'"/bin/sh" > user_setup.sh
echo "export PYTHIA8=$PWD" >> user_setup.sh
echo "export PYTHIA8DATA=\$PYTHIA8/share/Pythia8/xmldoc" >> user_setup.sh
echo "export LD_LIBRARY_PATH=\$PYTHIA8/lib:\$LD_LIBRARY_PATH" >> user_setup.sh
echo "export PYTHIA8_INCLUDE_DIR=\$PYTHIA8/include" >> user_setup.sh
echo "export PYTHIA8_LIBRARY=\$PYTHIA8/lib" >> user_setup.sh
source user_setup.sh
cd ..


echo "Installing Delphes 3.4.1"
#Check that the timing smearing is properly fixed tf_smeared
if [ ! -d "Delphes-3.4.1" ]; then
  # wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.4.1.tar.gz
  # tar xvzf Delphes-3.4.1.tar.gz
  # rm -rf Delphes-3.4.1.tar.gz
  # cd Delphes-3.4.1
  git clone https://github.com/ocerri/Delphes.git
  cd Delphes
fi
make -j8 HAS_PYTHIA8=true DelphesPythia8
cd ..

cd $start_dir
