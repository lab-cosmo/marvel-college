# marvel-college

## Docker

to build the Docker image:

docker build -t tutorial ./docker/

to open a notebook from the image:

docker run -it --rm  -v $HOME/:/home/jovyan -p 8888:8888 tutorial

(-v  set up a volume that links the ~/myhome:~/dokerhome )

## Notebooks

they are in the docker folder so that they can be copied to the image