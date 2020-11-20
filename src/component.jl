"""
This struct mean all atomic information other than coordinates.
This contains 5 fields like below.

- resname  : Union{String,  Nothing}
- resid    : Union{Int64,   Nothing}
- atomname : Union{String,  Nothing}
- atomid   : Union{Int64,   Nothing}
- mass     : Union{Float32, Nothing}
"""
mutable struct Attribute

    resname::Union{String, Nothing}
    resid::Union{Int64, Nothing}
    atomname::Union{String, Nothing}
    atomid::Union{Int64, Nothing}
    mass::Union{Float32, Nothing}

    function Attribute(;resname = nothing, resid = nothing,
                       atomname = nothing, atomid = nothing, mass = nothing)
        new(resname, resid, atomname, atomid, mass)
    end
end

"""
This struct have all information of corresponding atom.
This contains 2 fields like below.

- coordinate : Reference to [`Coordinate`](@ref) object.
- attribute  : Reference to [`Attribute`](@ref) object.
"""
mutable struct Atom
    coordinate::Coordinate{Float32}
    attribute::Attribute
end

function Atom(coordinate)
    Atom(coordinate, Attribute())
end

"""
This struct corresponding to MD trajectory.
This contains 4 fields like below.

- coordinates : AbstractArray{[`Coordinate`](@ref), 2}
- attributes  : Vector{[`Attribute`](@ref)}
- natom       : Number of atoms.
- nframe      : Number of frames in the trajectory.

The column of coordinates matrix means one snapshot.
The row of coordinates matrix means time series of one atom.

```julia
[ a d g j
  b e h k
  c f i l ]
```

In this case, one snapshot correspond to `[a b c]`.
"""
mutable struct Trajectory{RealT <: Real}
    coordinates::AbstractArray{Coordinate{RealT}, 2}
    attributes::Vector{Attribute}
    natom::Int64
    nframe::Int64

    function Trajectory(coordinates::AbstractArray{Coordinate{RealT}, 2},
                        attributes::AbstractArray{Attribute, 1} = [Attribute() for i=1:size(coordinates, 1)]) where RealT <: Real
        new{RealT}(coordinates, attributes, size(coordinates, 1), size(coordinates, 2))
    end
end

"""
This struct correspond to one snapshot of trajectory.
This contains 3 fields like below.

- coordinates : Vector{[`Coordinate`](@ref)}
- attribute   : Vector{[`Attribute`](@ref)}
- natom       : Number of atoms.
"""
mutable struct Frame{RealT <: Real}
    coordinates::Vector{Coordinate{RealT}}
    attributes::Vector{Attribute}
    natom::Int64

    function Frame(coordinates::Vector{Coordinate{RealT}},
                   attributes = [Attribute() for i=1:length(coordinates)]) where RealT <: Real
        new{RealT}(coordinates, attributes, length(coordinates))
    end
end
