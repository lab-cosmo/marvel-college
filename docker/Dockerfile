FROM jupyter/scipy-notebook:c7fb6660d096

COPY ./notebooks/ $HOME/notebooks/

RUN pip install ase==3.15
RUN pip install nglview==1.0
RUN jupyter-nbextension enable nglview --py --sys-prefix
