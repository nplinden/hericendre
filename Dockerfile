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
    zip \
    python3  \
    python3-dev


RUN git clone --recurse-submodules https://github.com/nplinden/hericendre.git && \
    cd hericendre && \
    mkdir build && \
    cd build && \
    cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE=Release .. && \
    make

CMD ["bash"]