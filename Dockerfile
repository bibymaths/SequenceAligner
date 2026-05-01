# ==========================================
# STAGE 1: Build the React Frontend
# ==========================================
FROM node:20-slim AS frontend-builder
WORKDIR /frontend

# Copy only package.json first to leverage Docker caching for npm install
COPY sequence_alignment_platform/frontend/package*.json ./
RUN npm install

# Copy the rest of the frontend source and build
COPY sequence_alignment_platform/frontend/ ./
RUN npm run build


# ==========================================
# STAGE 2: Build Backend & Final Image
# ==========================================
FROM python:3.12-slim
WORKDIR /app

# Install build dependencies (kept in final image due to C++ runtime requirements like openmpi)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    openmpi-bin \
    libopenmpi-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy the entire backend and C++ source
COPY . /app

# Build C++ tools
RUN mkdir -p /app/cpp_build && \
    cd /app/cpp_build && \
    cmake -DDOCKER_BUILD=ON .. && \
    make -j"$(nproc)"

# Install Python project
RUN pip install --upgrade pip && \
    pip install --no-cache-dir /app

# Fetch the compiled React static files from STAGE 1
RUN mkdir -p /app/sequence_alignment_platform/backend/static
COPY --from=frontend-builder /frontend/build/ /app/sequence_alignment_platform/backend/static/

EXPOSE 8000

CMD ["uvicorn", "sequence_alignment_platform.backend.main:app", "--host", "0.0.0.0", "--port", "8000"]