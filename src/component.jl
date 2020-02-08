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

mutable struct Atom

    coordinate::Vector{Float32} # [x, y, z]
    attribute::Attribute
end

function Atom(coordinate)
    Atom(coordinate, Attribute())
end

mutable struct Trajectory

    coordinates::Matrix{Vector{Float32}}
    # The column of coordinates matrix means one snapshot.
    # The row of coordinates matrix means time series of one atom.
    # [ a d g j
    #   b e h k
    #   c f i l ]
    # In this case, one snapshot correspond to [a b c]
    attributes::Vector{Attribute}
    nframe::Union{Int64, Nothing}

    function Trajectory(coordinates, attributes = [Attribute() for i=1:size(coordinates, 2)])
        new(coordinates, attributes, size(coordinates, 2))
    end
end
