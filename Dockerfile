# Copyright (C) 2025 Ladislav Hovan <ladislav.hovan@ncmbm.uio.no>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
# This library is free software: you can redistribute it and/or
# modify it under the terms of the GNU Public License as published
# by the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Public License along
# with this library. If not, see <https://www.gnu.org/licenses/>.

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
RUN python -m pip install --no-cache-dir netzoopy-sponge

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