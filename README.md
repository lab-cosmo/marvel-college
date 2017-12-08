# marvel-college

## Docker

to build the Docker image:

docker build -t tutorial ./docker/

to open a notebook from the image:

docker run -it --rm  -v ~/:~/ -p 8888:8888 tutorial

(-v  set up a volume that links the ~/myhome:~/dokerhome )