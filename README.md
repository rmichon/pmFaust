# Faust Physical Modeling ToolKit

This repository contains a set of tools to facilitate the implementation of physical models in the [Faust programming language](http://faust.grame.fr). This is an ongoing project funded in part the FEVER ANR project and carried out in partnership between [GRAME](http://grmae.fr), [CCRMA](https://ccrma.stanford.edu) (Stanford University) and [CRI](https://www.cri.mines-paristech.fr/) (Mines ParisTech).

## `pm.lib`

`pm.lib` is a physical modeling library for Faust. It contains various functions to conveniently connect bi-directional DSP blocks as well as higher level elements implementing different types of strings, tubes, solids, etc.

This is a work in progress and lots remains to be done.

## `IR2dsp.py`

`IR2dsp.py` is a simple python script taking an impulse response and converting it into a Faust modal physical model based on `pm.lib`. 

## `mesh2dsp.py`

`mesh2dsp.py` is a simple python script taking a 3D mesh and converting it into a Faust modal physical model based on `pm.lib`. 
