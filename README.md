# MolHandler

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://yutakasi634.github.io/MolHandler.jl/stable) -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://yutakasi634.github.io/MolHandler.jl/dev)
[![Build Status](https://travis-ci.com/yutakasi634/MolHandler.jl.svg?branch=master)](https://travis-ci.com/yutakasi634/MolHandler.jl)
[![Codecov](https://codecov.io/gh/yutakasi634/MolHandler.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yutakasi634/MolHandler.jl)
[![Coveralls](https://coveralls.io/repos/github/yutakasi634/MolHandler.jl/badge.svg?branch=master)](https://coveralls.io/github/yutakasi634/MolHandler.jl?branch=master)

## Concept
Simple and intuitive package for handling Molecular dynamics(MD) trajectory data.
Supported formats are `dcd` and `pdb`.

## How to install
This package assumes julia 1.2.
```julia
pkg> add https://github.com/yutakasi634/MolHandler.jl.git
julia> using MolHandler
```

## Example of use
```julia
julia> using MolHandler
julia> trj = read_dcd("trajectory.dcd")
julia> #trj = read_pdb("structure.pdb")
julia> trj.coordinates[:,1] # get first snapshot as Coordinate object array.
julia> trj.coordinates[1,:] # get first atom coordinate time series by Coordinate object array.
julia> frame      = get_frame(1, trj) # get first frame as Frame object.
julia> atom_array = get_atom(1, trj) # get first atom time series as Atom object array.
julia> atom       = get_atom(2, frame) # get second atom as Atom object.
```

## Document
Online document is https://yutakasi634.github.io/MolHandler.jl/dev.
