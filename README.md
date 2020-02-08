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
julia> #trj = readpdb("structure.pdb")
julia> trj.coordinates[:,1] # get first snapshot as Vector{Float32} array.
julia> trj.coordinates[1,:] # get first atom coordinate time series by Atom array.
julia> frame      = get_frame(1, trj) # get first frame as Frame object.
julia> atom_array = get_atom(1, trj) # get first atom time series as Atom array.
julia> atom       = get_atom(2, frame) # get second atom as Atom object.
```

### Trajectory struct
Trajectory struct have 2 fields like below.

    - coordinates : Matrix{Vector{Float32}}.
    - attributes  : Vector{Attribute}
    - nframe      : Number of frame in the trajectory.

### Frame struct
Frame struct have mean one snapshot of trajectory.
This struct have 2 fields like below.

    - coordinates : Vector{Vector{Float32}}
    - attribute   : Vector{Attribute}

### Atom struct
Atom struct have 2 fields like below.

    - coordinate : Vector{Float64}.
    - attribute  : Reference to Attribute class.

### Attribute struct
Attribute struct mean information of atom.
This struct have 5 fields like below.

    - resname  : Union{String,  Nothing}
    - resid    : Union{Int64,   Nothing}
    - atomname : Union{String,  Nothing}
    - atomid   : Union{Int64,   Nothing}
    - mass     : Union{Float32, Nothing}
