# Install file to create the environment of the EMStimTools library
# salome and gmsh are installed to /opt
#
# Author:
# Julius Zimmermann <julius.zimmermann@uni-rostock.de>

apt-get update
apt-get install -y software-properties-common

# Add Repository of 3rd Party Software, i.e. FEniCS
add-apt-repository ppa:fenics-packages/fenics

apt-get update
apt-get install -y \
    fenics \
    libgfortran3 \
    git \
    wget \
    bash-completion

wget http://security.ubuntu.com/ubuntu/pool/main/i/icu/libicu55_55.1-7ubuntu0.4_amd64.deb 2> /dev/null
wget http://security.ubuntu.com/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1.1_amd64.deb 2> /dev/null
dpkg -i /tmp/*.deb

wget https://bootstrap.pypa.io/get-pip.py 2> /dev/null
python3 get-pip.py
pip3 install --no-cache-dir setuptools

# Install SALOME
wget -O SALOME-8.3.0-UB16.04.tgz 'https://www.salome-platform.org/downloads/previous-versions/salome-v8.3.0/DownloadDistr?platform=UB16.04&version=8.3.0' 2> /dev/null
tar -xzf SALOME-8.3.0-UB16.04.tgz -C /opt

# Update path
export PATH="${PATH}:/opt/SALOME-8.3.0-UB16.04"

# Install GMSH
wget http://gmsh.info/bin/Linux/gmsh-4.0.2-Linux64-sdk.tgz 2> /dev/null
tar -xzf gmsh-4.0.2-Linux64-sdk.tgz -C /opt

# Update PYTHONPATH
export PYTHONPATH="/opt/gmsh-4.0.2-Linux64-sdk/lib"

# Clone and install EMSTimTools
git clone https://github.com/j-zimmermann/EMStimTools.git
cd EMStimTools || exit
python3 setup.py install
