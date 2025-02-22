
#!/bin/bash

# Compile the project
./build.sh

# Run test on pyramid.obj and output to pyramid.svg
./obj2cut pyramid.obj 10.0 10.0 > pyramid.svg

# Run test on head.obj and output to head.svg
./obj2cut head.obj 10.0 10.0 > head.svg