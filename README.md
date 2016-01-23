[![Build Status](https://travis-ci.org/jzuhone/formulas.svg?branch=master)](https://travis-ci.org/jzuhone/formulas)[![Coverage Status](https://coveralls.io/repos/github/jzuhone/formulas/badge.svg?branch=master)](https://coveralls.io/github/jzuhone/formulas?branch=master)

# formulas

Formulas is a package that serves as both a collection of common formulas in physics
and astrophysics and the machinery to create new formulas from scratch or from 
combinations of them. 

Formula objects are combinations of a SymPy-based mathematical formula with variables
and parameters, along with the ability to set the parameters to particular values. The
parameters may be set to float values or to unitful quantities, the latter provided by
the yt and/or AstroPy packages. 