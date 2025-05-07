ARG PYTHON_VERSION=3.12
FROM python:${PYTHON_VERSION}-slim

# Keeps Python from buffering stdout and stderr to avoid situations where
# the application crashes without emitting any logs due to buffering
ENV PYTHONUNBUFFERED=1

# Keeps Python from writing .pyc files
ENV PYTHONDONTWRITEBYTECODE=1

# Create a non-privileged user that the app will run under
# See https://docs.docker.com/go/dockerfile-user-best-practices/
ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --uid "${UID}" \
    --shell "/bin/bash" \
    appuser

# Install SPONGE
RUN python -m pip install netzoopy-sponge

# Clean up some space
RUN rm -rf /root/.cache/pip

# Create the working directory
WORKDIR /app
RUN chown appuser /app

# Switch to the non-privileged user to run the application
USER appuser

# Create an entry point
ENTRYPOINT ["netzoopy-sponge"]
CMD []

# Labels
LABEL org.opencontainers.image.source=https://github.com/kuijjerlab/sponge
LABEL org.opencontainers.image.description="Container image of SPONGE"
LABEL org.opencontainers.image.licenses=GPL-3.0-or-later