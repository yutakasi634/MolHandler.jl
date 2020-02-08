# MolHandler

[![Build Status](https://travis-ci.com/yutakasi634/MolHandler.svg?branch=master)](https://travis-ci.com/yutakasi634/MolHandler.jl)
[![Codecov](https://codecov.io/gh/yutakasi634/MolHandler.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yutakasi634/MolHandler.jl)
[![Coveralls](https://coveralls.io/repos/github/yutakasi634/MolHandler.jl/badge.svg?branch=master)](https://coveralls.io/github/yutakasi634/MolHandler.jl?branch=master)

## How to install
This package assumes julia 1.2.
```julia
pkg> add https://github.com/yutakasi634/MolHandler.jl.git
julia> using MolHandler
```

## How to use
```julia
julia> using MolHandler
julia> trj = readdcd("trajectory.dcd")
julia> trj.coordinates[:,1] # get first snapshot as Atom array.
julia> trj.coordinates[1,:] # get first atom coordinate time series by Atom array.
```

### Trajectory struct
Trajectory struct have 2 fields like below.

    - coordinates : Matrix{Atom}.
    - nframe      : Number of frame in the trajectory.

### Atom struct
Atom struct have 2 fields like below.

    - coordinate : Vector{Float64}.
    - attribute  : Reference to Attribute class. 

### Attribute struct
Attribute struct have 5 fields like below.

    - resname
    - resid
    - atomname
    - atomid
    - mass
