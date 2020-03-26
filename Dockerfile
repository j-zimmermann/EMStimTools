# Dockerfile to build the EMStimTools library

FROM ubuntu:18.04

WORKDIR /tmp

# needed since some other dependencies would otherwise try to read from the command line
ENV DEBIAN_FRONTEND=noninteractive

# Install requirements for software handling
RUN apt-get update && \
    apt-get install -y \
    software-properties-common

# Add Repository of 3rd Party Software
RUN add-apt-repository ppa:fenics-packages/fenics

# Install dependencies available via apt-get
# download needed libraries via wget
RUN apt-get update && \
    apt-get install -y \
    fenics \
    libgfortran3 \
    git \
    wget \
    bash-completion \
    libglu1-mesa-dev \
    python3-sphinx \
    latexmk \
    texlive-latex-extra && \
    apt-get clean

RUN wget http://security.ubuntu.com/ubuntu/pool/main/i/icu/libicu55_55.1-7ubuntu0.5_amd64.deb 2> /dev/null && \
    wget http://security.ubuntu.com/ubuntu/pool/main/libp/libpng/libpng12-0_1.2.54-1ubuntu1.1_amd64.deb 2> /dev/null && \
    dpkg -i /tmp/*.deb   && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install setuptools
RUN wget https://bootstrap.pypa.io/get-pip.py 2> /dev/null && \
    python3 get-pip.py && \
    pip3 install --no-cache-dir setuptools && \
    rm -rf /tmp/*

# Install SALOME
RUN wget -O SALOME-8.3.0-UB16.04.tgz 'https://www.salome-platform.org/downloads/previous-versions/salome-v8.3.0/DownloadDistr?platform=UB16.04&version=8.3.0' 2> /dev/null && \
    tar -xzf SALOME-8.3.0-UB16.04.tgz -C /opt && \
    rm -rf /tmp/*

# Update path
ENV PATH="${PATH}:/opt/SALOME-8.3.0-UB16.04"

# Install GMSH
RUN wget http://gmsh.info/bin/Linux/gmsh-4.0.2-Linux64-sdk.tgz 2> /dev/null && \
    tar -xzf gmsh-4.0.2-Linux64-sdk.tgz -C /opt && \
    rm -rf /tmp/*

# Update PYTHONPATH
ENV PYTHONPATH="/opt/gmsh-4.0.2-Linux64-sdk/lib"

# Clone and install EMSTimTools and python dependencies
RUN git clone https://github.com/j-zimmermann/EMStimTools.git 
RUN cd EMStimTools && \
    pip3 install .
