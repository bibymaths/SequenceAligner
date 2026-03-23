FROM python:3.11-slim AS base

# Install system dependencies needed for C++ compilation and analysis tools
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    gnuplot \
    imagemagick \
    jq \
    openmpi-bin \
    libopenmpi-dev \
    libdivsufsort-dev \
    git \
    nodejs \
    npm && \
    rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy full project so CMakeLists.txt at repo root is included
COPY . /app

# Build C++ tools from repo root
RUN mkdir -p /app/cpp_build && \
    cd /app/cpp_build && \
    cmake -DDOCKER_BUILD=ON .. && \
    make -j$(nproc)

# Copy Python backend
COPY sequence_alignment_platform/backend /app/backend

# Install Python dependencies
RUN pip install --no-cache-dir fastapi uvicorn[standard] python-multipart numpy

# Copy React frontend sources and build
COPY sequence_alignment_platform/frontend /app/frontend
RUN cd /app/frontend && \
    npm install && \
    npm run build

# Copy static files into backend static directory
RUN mkdir -p /app/backend/static && \
    cp -r /app/frontend/build/* /app/backend/static/

EXPOSE 8000

CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]