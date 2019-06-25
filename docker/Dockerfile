FROM jupyter/scipy-notebook:c7fb6660d096

# FROM felixmusil/appmode
# FROM appmode

# ARG NB_USER="jovyan"
# ARG NB_UID="1000"
# ARG NB_GID="100"



RUN pip install  ase==3.15 tabulate

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install -c bioconda nglview -y
RUN conda install -c conda-forge tqdm  -y
RUN conda install -c conda-forge  jupyter_contrib_nbextensions  -y
RUN conda install -c conda-forge   appmode -y

RUN jupyter-nbextension enable nglview --py --sys-prefix
RUN jupyter nbextension enable codefolding/main
RUN jupyter nbextension enable collapsible_headings/main
RUN jupyter nbextension enable toc2/main
RUN jupyter nbextension enable comment-uncomment/main
RUN jupyter nbextension enable keyboard_shortcut_editor/main
RUN jupyter nbextension enable codefolding/edit
RUN jupyter nbextension enable latex_envs/latex_envs
RUN jupyter nbextension enable equation-numbering/main
RUN jupyter serverextension enable --py --sys-prefix appmode
RUN jupyter nbextension     enable --py --sys-prefix appmode

USER root

RUN sudo apt-get update && \
    sudo apt-get install -y --no-install-recommends nodejs && \
    sudo apt-get install -y --no-install-recommends npm && \
    rm -rf /var/lib/apt/lists/*

USER $NB_UID

RUN cd /home/$NB_USER && \
    git clone https://github.com/osscar-org/widget-code-input.git && \
    cd widget-code-input && \
    pip install -e . && \
    jupyter nbextension install --py --symlink --sys-prefix widget_code_input && \
    jupyter nbextension enable --py --sys-prefix widget_code_input

RUN conda install -c conda-forge  matplotlib -y

COPY ./notebooks/ $HOME/notebooks/

USER root
# Launch Notebook server
EXPOSE 8888
CMD ["jupyter-notebook","--ip=127.0.0.1", "--allow-root", "--no-browser", "--NotebookApp.token=''"]
