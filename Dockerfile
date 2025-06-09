FROM python:3.10-slim-buster

# Install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc && \
    pip install --no-cache-dir pandas csvkit && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
