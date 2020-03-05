module MolHandler

export Coordinate, norm
export Atom, Attribute, Frame, Trajectory
export readdcd, readpdb
export get_frame, get_atom, clip_trajectory

# codes
include("coordinate.jl")
include("component.jl")
include("fileio.jl")
include("handling.jl")

end # module
