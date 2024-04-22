FROM ubuntu:jammy

ARG DEBIAN_FRONTEND=noninteractive 

RUN apt-get -y update && apt -y install \
    gcc git \
    wget vim curl \
    python3-pip cmake automake \
    build-essential \
    apt-utils flex bison mona \
    libeigen3-dev \
    coinor-clp coinor-libclp-dev

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

