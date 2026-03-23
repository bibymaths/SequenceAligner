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

# Copy full repository first so local Python packages and CMake root are available
COPY . /app

# Build C++ tools from repo root
RUN mkdir -p /app/cpp_build && \
    cd /app/cpp_build && \
    cmake -DDOCKER_BUILD=ON .. && \
    make -j"$(nproc)"

# Install Python dependencies
# Combined the requested packages with the existing ones, wrapping versions in quotes
RUN pip install --no-cache-dir \
    "numpy" \
    "pandas" \
    "matplotlib" \
    "mkdocs-material" \
    "pytest" \
    "mypy" \
    "ruff" \
    "mkdocs" \
    "mkdocstrings[python]" \
    "pystrings" \
    "tqdm>=4.67.3" \
    "scikit-learn>=1.8.0" \
    "umap-learn>=0.5.11" \
    "biopython>=1.86" \
    "fastapi>=0.135.1" \
    "uvicorn[standard]>=0.42.0" \
    "pydantic>=2.12.5" \
    "python-multipart>=0.0.22" \
    "websockets>=16.0" \
    "httpx>=0.28.1" && \
    if [ -d "/app/alignment_tool" ] || [ -f "/app/pyproject.toml" ] || [ -f "/app/setup.py" ]; then \
        pip install --no-cache-dir /app; \
    fi

# Copy Python backend
COPY sequence_alignment_platform/backend /app/backend

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