# Intrinsic Simplification
This is the minimum example to show how to use our intrinsic mesh simplification. To run this example, please compile it using the common cmake/make routine:
``` bash
cd 00_coarsening/
mkdir build
cd build
cmake ..
make -j8
```
Once compiled, one can run the example by typing
``` bash
./main
```
and you will see three examples
<img src="../assets/00.jpg" width="100%">
where the images (from left to right) are the original mesh, barycentric points of the original vertices on the coarsened mesh, and the coarsend mesh. We want to be crystally clear that the visualization of the coarsened mesh in the demo is INACCURATE. Speficially, the edge lengths in the visualization is different from the actual intrinsic edge lengths. The correct edge lengths require to extract the edge lengths matrix `l` from the code.

## Usage

You can simplify meshes by running the `./main` executable. By default, this simplifies the provided `bs_rest.obj` mesh down to 500 vertices. You can also specify a mesh and target coarseness as input by running
``` bash
./main /path/to/mesh.obj nVertices
```
The input mesh must be a manifold and connected obj file.
