FROM ubuntu:jammy

ARG DEBIAN_FRONTEND=noninteractive 

RUN apt-get -y update && apt -y install \
    gcc git \
    wget vim curl unzip swig \
    python3-pip cmake automake lsb-release \
    build-essential \
    apt-utils flex bison mona \
    libeigen3-dev 
    
WORKDIR /root/installation
RUN git clone https://github.com/google/or-tools
WORKDIR /root/installation/or-tools
RUN cmake -S . -B build -DBUILD_DEPS=ON 
#RUN cmake --build build --config Release --target all -j -v
RUN cmake --build build --config Release --target install -v


## Set environment variables
#ENV OR_TOOLS_VERSION="9.9"
#ENV OR_TOOLS_SUBVERSION="3963"
#ENV OR_TOOLS_URL="https://github.com/google/or-tools/releases/download/v${OR_TOOLS_VERSION}/or-tools_amd64_ubuntu-22.04_cpp_v${OR_TOOLS_VERSION}.${OR_TOOLS_SUBVERSION}.tar.gz"
#ENV OR_TOOLS_INSTALL_DIR="/or-tools"
#
## Download and extract OR-Tools
#RUN mkdir -p $OR_TOOLS_INSTALL_DIR && \
#    wget -O /tmp/or-tools.tar.gz $OR_TOOLS_URL && \
#    tar -xzf /tmp/or-tools.tar.gz -C $OR_TOOLS_INSTALL_DIR --strip-components=1 && \
#    rm -rf /tmp/or-tools.tar.gz
#
#RUN cd $OR_TOOLS_INSTALL_DIR && make test

## Install proxsuite
#WORKDIR /root/installations
#RUN git clone https://github.com/Simple-Robotics/proxsuite.git --recursive
#RUN mkdir build 
#WORKDIR /root/installations/proxsuite/build
#RUN cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DBUILD_WITH_VECTORIZATION_SUPPORT=OFF
#RUN make
#RUN make install


WORKDIR /root/berry-er

# Make the build directory if it doesnt exist, then remove it to make sure a no previous build files exist:
RUN mkdir -p build && rm -r build
RUN mkdir build 
#RUN cd build && cmake ./.. && make

