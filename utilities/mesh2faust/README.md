# mesh2faust

`mesh2faust` is command line tool to generate Faust (<http://faust.grame.fr>)
modal physical models from CAD files using finite element method (FEM). It is
based on [Vega FEM](http://run.usc.edu/vega/). `mesh2faust` should work both
on Linux and OSX.

## Build/Installation

`mesh2faust` relies on Vega FEM. A lightweight adapted version of this library
is part of this repository. In order to compile `mesh2faust`, some of Vega's
dependencies must be installed:

* Intel MKL library

To compile and install `mesh2faust`:

```
make
sudo make install
```

## Extra Help

Feel free to contact Romain Michon if you need extra help:
rmichon AT ccrma DOT stanford DOT edu.
