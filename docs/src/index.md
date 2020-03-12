# MolHandler.jl Documentation

This package is for handling MD trajectory data.
Supported format is `dcd` and `pdb`.

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
julia> atom_array = get_atom(1, trj) # get first atom time series as Atom array.
julia> atom       = get_atom(2, frame) # get second atom as Atom object.
```

## [Functions](@ref) list
```@index
Order = [:function]
```

## [Structs](@ref) list
```@index
Order = [:type]
```
