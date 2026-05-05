FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

# Install deps
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    python3 \
    python3-pip \
    python3-venv \
    autoconf \
    automake \
    libtool \
    libtool-bin \
    perl \
    flex \
    bison \
    ca-certificates \
    clang \
    libomp-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

# Invalidate old cache if needed
ADD https://api.github.com/repos/FUnTiDES-sim/FUnTiDES-TPL/commits/main /tmp/latest_commit

# set gcc 13 as build tool
ENV CC=clang
ENV CXX=clang++

# Clone TPL repo
RUN git clone --recurse-submodules https://github.com/FUnTiDES-sim/FUnTiDES-TPL.git tpl_source

# Install TPL
WORKDIR /workspace/tpl_source
RUN CC=clang CXX=clang++ ./install.sh --prefix=/opt/tpl --disable-cuda --enable-venv --jobs=$(nproc) --pykokkos-build-threads=4

ENV PATH="/opt/tpl/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/tpl/lib"
ENV CMAKE_PREFIX_PATH="/opt/tpl"

RUN rm -rf /workspace/tpl_source /tmp/latest_commit

WORKDIR /workspace
CMD ["/bin/bash"]
