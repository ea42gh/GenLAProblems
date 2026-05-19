FROM ea42gh/la-course:0.2

ENV JULIA_CONDAPKG_BACKEND=Null
ENV CONDAPKG_BACKEND=Null
ENV JULIA_PYTHONCALL_EXE=/usr/local/bin/python3
ENV JULIA_PROJECT=/home/jovyan/GenLAProblems
ENV JULIACALL_PROJECT=/home/jovyan/GenLAProblems
ENV PYTHON_JULIACALL_EXE=/usr/local/julia/bin/julia
ENV PYTHON_JULIACALL_PROJECT=/home/jovyan/GenLAProblems
ENV PYJULIAPKG_PROJECT=/home/jovyan/GenLAProblems
ENV JULIAPKG_PROJECT=/home/jovyan/GenLAProblems
ENV JULIA_PKG_PRECOMPILE_AUTO=0
ENV JULIA_NUM_THREADS=1
ENV JULIA_NUM_PRECOMPILE_TASKS=1
ENV OPENBLAS_NUM_THREADS=1
ENV MKL_NUM_THREADS=1
ENV OMP_NUM_THREADS=1
ENV TERM=xterm-256color
ENV INPUTRC=/home/jovyan/.inputrc
ENV XDG_DATA_HOME=/home/jovyan/.local/share
ENV JUPYTER_DATA_DIR=/home/jovyan/.local/share/jupyter

USER root
WORKDIR /tmp

RUN mkdir -p /home/jovyan/GenLAProblems \
 && chown -R jovyan:jovyan /home/jovyan/GenLAProblems

RUN /usr/local/bin/python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel \
 && /usr/local/bin/python3 -m pip install --no-cache-dir \
    sympy \
    "matrixlayout[render] @ git+https://github.com/ea42gh/matrixlayout.git@e65aa0ab6e438b8e4dc961a0dae5f4360c480cd6" \
    "LAFigureSpecs @ git+https://github.com/ea42gh/LAFigureSpecs.git@84bec61d56f07fffa4ca7e9c0476c12e2847e274"

USER jovyan
WORKDIR /home/jovyan/GenLAProblems

COPY . /home/jovyan/GenLAProblems

USER root
RUN chown -R jovyan:jovyan /home/jovyan/GenLAProblems

RUN mkdir -p /home/jovyan/work \
 && chown -R jovyan:jovyan /home/jovyan/work \
 && chmod -R u+rwX /home/jovyan/work

COPY --chown=jovyan:jovyan GenProblems.ipynb /home/jovyan/work/GenProblems.ipynb
COPY --chown=jovyan:jovyan KernelProbe.ipynb /home/jovyan/work/KernelProbe.ipynb

USER jovyan

RUN bash binder/postBuild

EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--ServerApp.token=", "--ServerApp.password=", "--ServerApp.root_dir=/home/jovyan/work", "--notebook-dir=/home/jovyan/work"]
