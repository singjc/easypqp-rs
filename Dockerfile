# Multi-stage build with CUDA-enabled base images so the final image can use
# the crate's `cuda` feature. The builder installs Rust via rustup and builds
# the `easypqp-insilico` binary with the `cuda` feature enabled.

# Stage 1: Builder with CUDA dev image
FROM nvidia/cuda:12.2.2-devel-ubuntu22.04 AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    libssl-dev \
    pkg-config \
    git \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Rust toolchain (rustup)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# CUDA environment variables for build
ENV CUDA_HOME=/usr/local/cuda
ENV PATH=${CUDA_HOME}/bin:${PATH}
ENV LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}

# Work directory and copy sources
WORKDIR /build
COPY Cargo.toml ./
COPY easypqp-core ./easypqp-core
COPY easypqp-cli ./easypqp-cli
COPY easypqp-py ./easypqp-py

# Build release binary with CUDA feature enabled
RUN cargo build --release -p easypqp-cli --features cuda

# Stage 2: Runtime with CUDA runtime image
FROM nvidia/cuda:12.2.2-runtime-ubuntu22.04 AS runtime

# Install minimal runtime dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    libstdc++6 \
    libgomp1 \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy binary from builder
COPY --from=builder /build/target/release/easypqp-insilico /usr/local/bin/

# CUDA env for runtime
ENV CUDA_HOME=/usr/local/cuda
ENV PATH=${CUDA_HOME}/bin:${PATH}
ENV LD_LIBRARY_PATH=${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}

WORKDIR /data

ENTRYPOINT ["easypqp-insilico"]
CMD ["--help"]
