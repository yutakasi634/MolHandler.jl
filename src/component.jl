struct Atom

    coordinate::Vector{Float64} # [x, y, z]
    resname::String
    resid::Int64
    atomname::String
    atomid::Int64
    mass::Float64

    Atom(coordinate) = new(coordinate)
end

struct Trajectory

    coordinates::Matrix{Atom}
    # The column of coordinates matrix means one snapshot.
    # The row of coordinates matrix means time series of one atom.
    # [ a d g j
    #   b e h k
    #   c f i l ]
    # In this case, one snapshot correspond to [a b c]
    nframe::Int64

    Trajectory(coordinates) = new(coordinates, size(coordinates, 2))
end
