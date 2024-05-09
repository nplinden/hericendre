FROM --platform=linux/x86-64 ubuntu:22.04

WORKDIR /root

RUN apt update -y && \
    apt install -y \ 
    wget \
    git \
    zsh \
    cmake \
    gcc \
    g++ \
    libeigen3-dev \
    libpugixml-dev \
    libyaml-cpp-dev \
    libyaml-cpp0.7 \
    zip \
    python3 

RUN wget https://github.com/fmtlib/fmt/releases/download/10.2.1/fmt-10.2.1.zip && \
    unzip fmt-10.2.1.zip && \
    cd fmt-10.2.1 && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install

RUN git clone --recurse-submodules https://github.com/nplinden/hericendre.git && \
    cd hericendre && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release .. && \
    make

CMD ["bash"]