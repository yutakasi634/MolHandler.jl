# MolHandler

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://yutakasi634.github.io/MolHandler.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://yutakasi634.github.io/MolHandler.jl/dev)
[![Build Status](https://travis-ci.com/yutakasi634/MolHandler.svg?branch=master)](https://travis-ci.com/yutakasi634/MolHandler.jl)
[![Codecov](https://codecov.io/gh/yutakasi634/MolHandler.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yutakasi634/MolHandler.jl)
[![Coveralls](https://coveralls.io/repos/github/yutakasi634/MolHandler.jl/badge.svg?branch=master)](https://coveralls.io/github/yutakasi634/MolHandler.jl?branch=master)

## How to install
This package assumes julia 1.2.
```julia
pkg> add https://github.com/yutakasi634/MolHandler.jl.git
julia> using MolHandler
```

## Example of use
```julia
julia> using MolHandler
julia> trj = readdcd("trajectory.dcd")
julia> #trj = readpdb("structure.pdb")
julia> trj.coordinates[:,1] # get first snapshot as Vector{Float32} array.
julia> trj.coordinates[1,:] # get first atom coordinate time series by Atom array.
julia> frame      = get_frame(1, trj) # get first frame as Frame object.
julia> atom_array = get_atom(1, trj) # get first atom time series as Atom array.
julia> atom       = get_atom(2, frame) # get second atom as Atom object.
```

## Document
Online document is https://yutakasi634.github.io/MolHandler.jl/dev.
