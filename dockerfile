# Use official micromamba base image (recommended by recent nf-core guidelines)
FROM mambaorg/micromamba:1.5.3

# Metadata labels (Best Practice)
LABEL image.author="[Your Name]" \
      image.description="Docker image containing all requirements for nf-core/manifold pipeline" \
      image.source="https://github.com/nf-core/manifold" \
      image.licenses="MIT"

# Copy the environment file into the container
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/env.yaml

# Install dependencies
# Note: "clean -a" drastically reduces the final image size
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean -a -y

# Add conda environment to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Default command
CMD ["/bin/bash"]
