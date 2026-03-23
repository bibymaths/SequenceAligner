FROM python:3.12-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    openmpi-bin \
    libopenmpi-dev \
    libdivsufsort-dev \
    git \
    nodejs \
    npm && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . /app

# Build C++ tools
RUN mkdir -p /app/cpp_build && \
    cd /app/cpp_build && \
    cmake -DDOCKER_BUILD=ON .. && \
    make -j"$(nproc)"

# Install Python project
RUN python -m pip install --upgrade pip setuptools wheel && \
    pip install --no-cache-dir /app

# Build React frontend
RUN cd /app/sequence_alignment_platform/frontend && \
    npm install && \
    npm run build

# Copy build → backend/static
RUN mkdir -p /app/sequence_alignment_platform/backend/static && \
    cp -r /app/sequence_alignment_platform/frontend/build/. \
          /app/sequence_alignment_platform/backend/static/

EXPOSE 8000

CMD ["uvicorn", "sequence_alignment_platform.backend.main:app", "--host", "0.0.0.0", "--port", "8000"]