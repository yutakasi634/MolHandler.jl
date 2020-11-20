# MolHandler

<!-- [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://yutakasi634.github.io/MolHandler.jl/stable) -->
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://yutakasi634.github.io/MolHandler.jl/dev)
[![Build Status](https://travis-ci.com/yutakasi634/MolHandler.jl.svg?branch=master)](https://travis-ci.com/yutakasi634/MolHandler.jl)
[![Codecov](https://codecov.io/gh/yutakasi634/MolHandler.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/yutakasi634/MolHandler.jl)
[![Coveralls](https://coveralls.io/repos/github/yutakasi634/MolHandler.jl/badge.svg?branch=master)](https://coveralls.io/github/yutakasi634/MolHandler.jl?branch=master)

## Concept
Simple and intuitive package for handling Molecular dynamics trajectory data.
Supported formats are `xyz`, `pdb` and `dcd`.

## How to install
Julia version 1.2 is supported.
```julia
pkg> add https://github.com/yutakasi634/MolHandler.jl.git
julia> using MolHandler
```

## Example of use
### Handling dcd file
```julia
julia> using MolHandler
julia> trj = read_dcd("trajectory.dcd")
julia> trj.coordinates[:,1] # get first snapshot as Coordinate object array.
julia> trj.coordinates[1,:] # get first atom coordinate time series by Coordinate object array.
julia> frame            = get_frame(1, trj)  # get first frame as Frame object.
julia> atom_time_series = get_atom(1, trj)   # get first atom time series as Atom object array.
julia> atom             = get_atom(2, frame) # get second atom as Atom object.
julia> atom.coordinate.x # return x element of atom  coordinate.
julia> com = center_of_mass(trj, geometric = ture)
julia> trj.coordinates .-= com[1] # move this trajectory to origin by first frame center of mass.
julia> write_dcd("com_fixed_trj.dcd", trj) # output this protein as dcd file.
```

### Handling pdb file
```julia
julia> using MolHandler
julia> trj = read_pdb("1aki.pdb", model = :AA) # read protein as all atom model.
julia> com = center_of_mass(trj)
julia> rg  = radius_of_gyration(trj)
julia> trj.coordinates .-= com # move this protein to origin.
julia> write_pdb("com_fixed_1aki.pdb", trj) # output this protein as pdb file.
```

### Support Multi-Threading
```julia
julia> using Molhandler
julia> trj = read_dcd("trajectory.dcd")
julia> contact_map = contact_probability_matrix_parallel(10.0, trj) # get contact-map of which contact threshold is 10.0.
```
You have to set the environment variable `JULIA_NUM_THREADS` before excution.

# Document
Online document is https://yutakasi634.github.io/MolHandler.jl/dev.
