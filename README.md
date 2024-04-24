This is a compilation of algorithms that evaluate Bézier curves.

# Requirements:
- mpfr
- eigen

# Compilation
Download the codes and extract. In the main folder, follow the instructions.

On Windows
~~~~
mkdir build
cd build
cmake .. -G Ninja
cmake --build .
~~~~
On Ubuntu
~~~~
mkdir build
cd build
cmake ..
cmake --build .
~~~~

# Manual
On Windows
```
.\compare.exe              get the runtime on a random curve for 1000 evaluation points
.\compare.exe -rd          get the runtime with respect to the degrees
.\compare.exe -rs          get the runtime with respect to the number of evaluation points
.\compare.exe -e [t]       get the relative errors at a give t
.\compare.exe -e           get the relative errors for an entire curve (t in [0,1])
```
On Ubuntu
```
./compare                  get the runtime on a random curve for 1000 evaluation points
./compare -rd              get the runtime with respect to the degrees
./compare -rs              get the runtime with respect to the number of evaluation points
./compare -e [t]           get the relative errors at a give t
./compare -e               get the relative errors for an entire curve (t in [0,1])
```
