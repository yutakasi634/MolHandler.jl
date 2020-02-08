module MolHandler

export Atom, Attribute, Frame, Trajectory
export readdcd, readpdb
export get_frame, get_atom

# codes
include("component.jl")
include("fileio.jl")
include("handling.jl")

end # module
