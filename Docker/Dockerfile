# Use a common Python base image suitable for both scripts
FROM python:3.12-slim

# Ensure apt is up-to-date and install system packages required for science
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get update && \
    apt-get install -y \
    build-essential \
    libopenblas-dev \
    python3-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set up Python environment
WORKDIR /app

# Copy requirements into image
COPY requirements.txt requirements.txt

# Install dependencies via pip
RUN pip install --upgrade pip && \
    pip install -r requirements.txt

# Optional: Copy your code (if placing scripts here directly)
# COPY ipynb/CNN_Local_noAutoencoders_3x3.py ipynb/CNN_Local_noAutoencoders_3x3.py
# COPY ipynb/CNN_Federated_noAutoencoder_3x3.py ipynb/CNN_Federated_noAutoencoder_3x3.py

# Entrypoint — override as needed, for interactive use just use: docker run -it ...
CMD ["/bin/bash"]