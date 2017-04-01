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

## License

Copyright 2017, Romain Michon and Sara R. Martin

Project partly funded by Research Council of Norway et NTNU (Norwegian Technical University of Trondheim).

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
