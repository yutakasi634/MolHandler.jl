var documenterSearchIndex = {"docs":
[{"location":"functions/#Functions-1","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/#Index-1","page":"Functions","title":"Index","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"Order = [:function]","category":"page"},{"location":"functions/#File-IO-1","page":"Functions","title":"File IO","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"readdcd\nreadpdb","category":"page"},{"location":"functions/#MolHandler.readdcd","page":"Functions","title":"MolHandler.readdcd","text":"readdcd(filename::String)::Trajectory\n\nReturn Trajectory object which filled coordinates, nframe, natom fields.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.readpdb","page":"Functions","title":"MolHandler.readpdb","text":"readpdb(filename::String)::Trajectory\n\nReturn Trajectory object which filled all field.\n\n\n\n\n\n","category":"function"},{"location":"functions/#Trajectory-handling-1","page":"Functions","title":"Trajectory handling","text":"","category":"section"},{"location":"functions/#","page":"Functions","title":"Functions","text":"get_frame\nget_atom","category":"page"},{"location":"functions/#MolHandler.get_frame","page":"Functions","title":"MolHandler.get_frame","text":"get_frame(frame_idx::Int64, trajectory::Trajectory)::Frame\n\nReturn Frame object which correspond to framd_idx frame from trajectory.\n\n\n\n\n\n","category":"function"},{"location":"functions/#MolHandler.get_atom","page":"Functions","title":"MolHandler.get_atom","text":"get_atom(query::Int64, trajectory::Trajectory)::Vector{Atom}\n\nReturn Vector of queryth Atom object from trajectory.\n\n\n\n\n\nget_atom(query::Int64, frame::Frame)::Atom\n\nReturn queryth Atom object from frame.\n\n\n\n\n\n","category":"function"},{"location":"structs/#Structs-1","page":"Structs","title":"Structs","text":"","category":"section"},{"location":"structs/#Index-1","page":"Structs","title":"Index","text":"","category":"section"},{"location":"structs/#","page":"Structs","title":"Structs","text":"Order = [:type]","category":"page"},{"location":"structs/#","page":"Structs","title":"Structs","text":"Trajectory\nCoordinate\nAttribute\nFrame\nAtom","category":"page"},{"location":"structs/#MolHandler.Trajectory","page":"Structs","title":"MolHandler.Trajectory","text":"This struct corresponding to MD trajectory. This contains 4 fields like below.\n\ncoordinates : Matrix{Coordinate}\nattributes  : Vector{Attribute}\nnatom       : Number of atoms.\nnframe      : Number of frames in the trajectory.\n\nThe column of coordinates matrix means one snapshot. The row of coordinates matrix means time series of one atom.\n\n[ a d g j\n  b e h k\n  c f i l ]\n\nIn this case, one snapshot correspond to [a b c].\n\n\n\n\n\n","category":"type"},{"location":"structs/#MolHandler.Coordinate","page":"Structs","title":"MolHandler.Coordinate","text":"This correspond to coordinate of atom. This contain 3 fields like below.\n\nx <: Real\ny <: Real\nz <: Real\n\nSome operators and functions were overloaded to make Coordinate objects broadcastable.\n\n\n\n\n\n","category":"type"},{"location":"structs/#MolHandler.Attribute","page":"Structs","title":"MolHandler.Attribute","text":"This struct mean all atomic information other than coordinates. This contains 5 fields like below.\n\nresname  : Union{String,  Nothing}\nresid    : Union{Int64,   Nothing}\natomname : Union{String,  Nothing}\natomid   : Union{Int64,   Nothing}\nmass     : Union{Float32, Nothing}\n\n\n\n\n\n","category":"type"},{"location":"structs/#MolHandler.Frame","page":"Structs","title":"MolHandler.Frame","text":"This struct correspond to one snapshot of trajectory. This contains 3 fields like below.\n\ncoordinates : Vector{Coordinate}\nattribute   : Vector{Attribute}\nnatom       : Number of atoms.\n\n\n\n\n\n","category":"type"},{"location":"structs/#MolHandler.Atom","page":"Structs","title":"MolHandler.Atom","text":"This struct have all information of corresponding atom. This contains 2 fields like below.\n\ncoordinate : Reference to Coordinate object.\nattribute  : Reference to Attribute object.\n\n\n\n\n\n","category":"type"},{"location":"#MolHandler.jl-Documentation-1","page":"Top","title":"MolHandler.jl Documentation","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"This package is for handling MD trajectory data. Supported format is dcd and pdb.","category":"page"},{"location":"#How-to-install-1","page":"Top","title":"How to install","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"This package assumes julia 1.2.","category":"page"},{"location":"#","page":"Top","title":"Top","text":"pkg> add https://github.com/yutakasi634/MolHandler.jl.git\njulia> using MolHandler","category":"page"},{"location":"#Example-of-use-1","page":"Top","title":"Example of use","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"julia> using MolHandler\njulia> trj = readdcd(\"trajectory.dcd\")\njulia> #trj = readpdb(\"structure.pdb\")\njulia> trj.coordinates[:,1] # get first snapshot as Coordinate object array.\njulia> trj.coordinates[1,:] # get first atom coordinate time series by Coordinate object array.\njulia> frame      = get_frame(1, trj) # get first frame as Frame object.\njulia> atom_array = get_atom(1, trj) # get first atom time series as Atom array.\njulia> atom       = get_atom(2, frame) # get second atom as Atom object.","category":"page"},{"location":"#[Functions](@ref)-list-1","page":"Top","title":"Functions list","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"Order = [:function]","category":"page"},{"location":"#[Structs](@ref)-list-1","page":"Top","title":"Structs list","text":"","category":"section"},{"location":"#","page":"Top","title":"Top","text":"Order = [:type]","category":"page"}]
}
