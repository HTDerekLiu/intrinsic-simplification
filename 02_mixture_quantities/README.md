# Mix-and-match Coarsening Quantities 
Our coarsening supports different "masses" for decimation. In this demo, we show how can one mix-and-match the curvature quantity (encourging curvature aware decimation) with the area quantity (encouraging uniform simplification) to achieve different decimation behaviors. 

To run this example, please compile it using the common cmake/make routine:
``` bash
cd 02_mixture_quantities/
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
<img src="../assets/02.jpg" width="100%">
From left to right, the results are (Gaussian) curvature only simplification, mixture of curvature and area simplification, area only simplification. Again, the visualization does not reflect the actualy intrinsic mesh. One should use the edge lengths matrix `l` in the code to extract the actual intrinsic geometry. 

## Usage

You can simplify meshes by running the `./main` executable. By default, this simplifies the provided `capsule.obj` mesh down to 500 vertices. You can also specify a mesh and target coarseness as input by running
``` bash
./main /path/to/mesh.obj nVertices
```
The input mesh must be a manifold and connected obj file.
