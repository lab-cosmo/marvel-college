FROM jupyter/scipy-notebook:c7fb6660d096

COPY ./notebooks/ $HOME/notebooks/

RUN pip install  ase==3.15

RUN conda config --add channels conda-forge
RUN conda install -c bioconda nglview -y
RUN conda install -c conda-forge tqdm -y

RUN conda install -c conda-forge jupyter_contrib_nbextensions -y
RUN jupyter-nbextension enable nglview --py --sys-prefix
RUN jupyter nbextension enable codefolding/main
RUN jupyter nbextension enable collapsible_headings/main
RUN jupyter nbextension enable toc2/main
RUN jupyter nbextension enable comment-uncomment/main
RUN jupyter nbextension enable keyboard_shortcut_editor/main
RUN jupyter nbextension enable codefolding/edit
RUN jupyter nbextension enable latex_envs/latex_envs
RUN jupyter nbextension enable equation-numbering/main
