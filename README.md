# marvel-college

## Docker

to build the Docker image:

docker build -t tutorial ./docker/

to open a notebook from the image:

docker run -it --rm -p 8888:8888 -v /local/git/marvel-college/MD/:/home/jovyan/notebooks/ tutorial

/local/git/marvel-college/MD/ should be changed to the path where the notebooks are located.

## Notebooks

they are in the docker folder so that they can be copied to the image