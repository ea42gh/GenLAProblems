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
ENV TERM=xterm-256color
ENV INPUTRC=/home/jovyan/.inputrc
ENV XDG_DATA_HOME=/home/jovyan/.local/share
ENV JUPYTER_DATA_DIR=/home/jovyan/.local/share/jupyter

USER root
WORKDIR /tmp

RUN mkdir -p /home/jovyan/GenLAProblems \
 && chown -R jovyan:jovyan /home/jovyan/GenLAProblems

RUN /usr/local/bin/python3 -m pip install --no-cache-dir --upgrade pip setuptools wheel sympy

USER jovyan
WORKDIR /home/jovyan/GenLAProblems

COPY . /home/jovyan/GenLAProblems

USER root
RUN chown -R jovyan:jovyan /home/jovyan/GenLAProblems

RUN mkdir -p /home/jovyan/work \
 && chown -R jovyan:jovyan /home/jovyan/work

COPY --chown=jovyan:jovyan GenProblems.ipynb /home/jovyan/work/GenProblems.ipynb

USER jovyan

RUN julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); using GenLAProblems, PythonCall; GenLAProblems.import_sympy()'
RUN bash binder/postBuild

EXPOSE 8888
CMD ["jupyter", "lab", "--ip=0.0.0.0", "--no-browser", "--ServerApp.token=", "--ServerApp.password=", "--ServerApp.root_dir=/home/jovyan/work", "--notebook-dir=/home/jovyan/work"]
