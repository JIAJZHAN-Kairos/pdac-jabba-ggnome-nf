FROM --platform=linux/amd64 mambaorg/micromamba:2.3.3@sha256:800e7ade3ffe29c9a9ac2026163131495f8197c3852e572c5835beb4e8a33cd6 AS mamba
FROM --platform=linux/amd64 mskilab/jabba@sha256:be8e391cf9e13c2870e4578fa80433315cad2f8c23dee5c55a1cdce660f54d47
USER root
ADD runtime_pinned.tar.gz /
ADD cplex_include_lib_original.tar.gz /
RUN test -s /g/data/gx8/projects/Jiajun_intogen/software/jz9027_CPLEX_Studio2220/cplex/include/ilcplex/cplex.h \
    && test -s /g/data/gx8/projects/Jiajun_intogen/software/jz9027_CPLEX_Studio2220/cplex/lib/x86-64_linux/static_pic/libcplex.a
COPY --from=mamba /bin/micromamba /usr/local/bin/micromamba
RUN micromamba create -y -p /opt/python -c conda-forge python=3.12.8 numpy=1.26.4 pandas=2.2.2 matplotlib-base=3.9.2 pip=24.3.1 \
    && /opt/python/bin/python -m pip install --no-cache-dir awscli==1.36.40 \
    && micromamba clean --all --yes
ENV PATH=/opt/python/bin:$PATH
ENV CPLEX_DIR=/g/data/gx8/projects/Jiajun_intogen/software/jz9027_CPLEX_Studio2220 \
    R_LIBS_USER=/g/data/gx8/projects/Jiajun_intogen/software/jz9027_JaBbA_R_20261027 \
    GUTILS_COMPAT_LIB=/g/data/gx8/projects/Jiajun_intogen/software/jz9027_gUtils_R40_compat_20261028 \
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/tmp/matplotlib

ENTRYPOINT []
