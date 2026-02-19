# Multi-stage build for easypqp-insilico binary
FROM rust:1.88-slim AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /build

# Copy workspace files
COPY Cargo.toml ./
COPY easypqp-core ./easypqp-core
COPY easypqp-cli ./easypqp-cli
COPY easypqp-py ./easypqp-py

# Build release binary with default features (includes parquet)
RUN cargo build --release -p easypqp-cli

# Runtime stage - minimal image
FROM debian:bookworm-slim

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    ca-certificates \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy binary from builder
COPY --from=builder /build/target/release/easypqp-insilico /usr/local/bin/

# Set working directory
WORKDIR /data

# Run the binary
ENTRYPOINT ["easypqp-insilico"]
CMD ["--help"]
