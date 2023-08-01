# Intrinsic Prolongation
This is an example to show how to construct a prolongation operator (a.k.a. interpolation) from our intrinsically simplified mesh. To run this example, please compile it using the common cmake/make routine:
```
cd 01_prolongation_scalar/
mkdir build
cd build
cmake ..
make -j8
```
Once compiled, one can run the example by typing
```
./main
```
and you will see two examples
<img src="../assets/01.jpg" width="100%">
where the first image is the test function (shown as colors) defined on the  simplified model and the second image is the prolonged counterpart of the test function on the original model. Again, we want to emphasize that the visualization of the simplified model is INACCURATE, meaning the straight edge length you see in the demo does not have the same length as the actual intrinsic edge length.