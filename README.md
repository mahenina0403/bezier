Requirements:

`mpfr`

How to compile?
> ```
mkdir build
cd build
cmake ..
make
>```

How to use it?
To just use the default benchmarking

`./compare -t`

To compare with a specific number of samples

`./compare -t [SAMPLES]`

example: `./compare -t 500`

To compare the relative error to the output of the de Casteljau algorithm

`./compare -p`

One can use different parameter by adding the parameter value

`./compare -p [parameter_value]`

example: `./compare -t 0.3`
