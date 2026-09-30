FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        build-essential \
        libgomp1 \
        time \
        coreutils \
        procps \
        python3-minimal && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /work

CMD ["/bin/bash"]