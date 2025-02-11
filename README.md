# obj2cut

obj2cut is a command-line tool that converts an OBJ file into an SVG file. It projects and packs the triangles from the 3D model into a 2D layout and draws both outer and inner triangles with labeled edges.

## Features
- Reads OBJ files with vertices and faces.
- Projects and packs triangles into a 2D space.
- Computes an inner triangle with a configurable offset.
- Generates an SVG with labeled edges.

## Usage
```bash
./obj2cut path/to/file.obj [scale (default:10.0)] [innerOffset (default:10.0)]
```

Example:
```bash
./obj2cut pyramid.obj 10.0 10.0 > output.svg
```

## Build Instructions
Compile using the following command:
```bash
g++ -std=c++11 -o obj2cut obj2cut.cc
```

## Notes
- The tool outputs the SVG to standard output.
  