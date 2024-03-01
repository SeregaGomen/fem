FROM rust
#FROM ubuntu

WORKDIR /
COPY . /.

# Update and upgrade repo
RUN DEBIAN_FRONTEND=noninteractive apt-get -y update

# Install tools we might need
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install curl

# Install Rust
RUN curl -y --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

#Install libs
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install liblapacke-dev libmumps-seq-dev libopenblas-dev libsuitesparse-dev

RUN apt -y install cargo
#RUN cargo build

RUN cargo build --release --bin example1
CMD ["./target/release/example1"]

#https://stackoverflow.com/questions/41092587/passing-a-file-as-an-argument-to-a-docker-container
