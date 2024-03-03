FROM rust:1.67

WORKDIR /
COPY . /.

# Update and upgrade repo
RUN DEBIAN_FRONTEND=noninteractive apt-get -y update

#Install libs
RUN DEBIAN_FRONTEND=noninteractive apt-get -y install liblapacke-dev libmumps-seq-dev libopenblas-dev libsuitesparse-dev

# Build project
RUN cargo build --release --bin example1


CMD ["./target/release/example1"]
