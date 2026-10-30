ARG MINIFORGE_VERSION=26.1.1-2
ARG UBUNTU_VERSION=24.04
ARG CONDA_ENV_PATH=/opt/conda/envs/proteovault

FROM condaforge/miniforge3:${MINIFORGE_VERSION} AS builder

ARG CONDA_ENV_PATH

RUN mamba create -qy -p "${CONDA_ENV_PATH}" \
        -c conda-forge \
        python=3.12 \
        pip && \
    mamba clean -afy

WORKDIR /build
COPY . .
RUN "${CONDA_ENV_PATH}/bin/python" -m pip install --no-cache-dir . \
    && "${CONDA_ENV_PATH}/bin/proteovault" --help > /dev/null

FROM ubuntu:${UBUNTU_VERSION} AS final

ARG CONDA_ENV_PATH

COPY --from=builder "${CONDA_ENV_PATH}" "${CONDA_ENV_PATH}"

ENV CONDA_ENV_PATH="${CONDA_ENV_PATH}" \
    PATH="${CONDA_ENV_PATH}/bin:${PATH}"

RUN groupadd -g 500001 bldocker \
    && useradd -r -u 500001 -g bldocker bldocker

USER bldocker
WORKDIR /work

LABEL maintainer="Chenghao Zhu <czhu@sbpdiscovery.org>" \
      org.opencontainers.image.source="https://github.com/TheBoutrosLab/package-ProteoVault" \
      org.opencontainers.image.title="ProteoVault" \
      org.opencontainers.image.description="Portable proteomics analysis results in DuckDB"
