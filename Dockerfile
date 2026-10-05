# Container for the whole pipeline.
#
#   docker build -t alnfail .
#   docker run --rm -it -v "$PWD":/work alnfail \
#       snakemake --cores 4 --configfile config/test.yaml
#
# For a locked build, export environment.lock.yml (see environment.yml) and
# change the COPY line below to use it.
FROM mambaorg/micromamba:1.5.10

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# make the environment's tools available to plain `docker run ... <command>`
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH
WORKDIR /work
