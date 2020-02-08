module MolHandler

export Atom, Attribute, Frame, Trajectory
export readdcd
export get_frame, get_atom, display_frame

# codes
include("component.jl")
include("fileio.jl")
include("handling.jl")

end # module
