import Dates
import Printf

function read_dcd_meta_info(io::IOStream)::DCDMetaInfo
    unitcell_flag = false # for charmm format

    seekend(io)
    file_size = position(io)
    seekstart(io)

    # read header first block
    skip(io, 4) # skip block size part
    header_sig = Array{Char, 1}(undef, 4)
    for i in 1:4
        header_sig[i] = read(io, Char)
    end
    total_frame_in_header = read(io, Int32)
    first_step = read(io, Int32)
    nstep_save = read(io, Int32)
    total_step = read(io, Int32)
    total_unit = read(io, Int32)
    header_null_4 = Array{Int32, 1}(undef, 4)
    read!(io, header_null_4)
    time_step  = read(io, Float32)
    header_null_9 = Array{Int32, 1}(undef, 9)
    read!(io, header_null_9)
    if header_null_9[1] != 0
        unitcell_flag = true
    end
    version   = read(io, Int32)
    skip(io, 4) # skip block size part

    # read header second block
    skip(io, 4)
    number_of_lines    = read(io, Int32)
    title = []
    for i in 1:number_of_lines
        line = Array{Char, 1}(undef, 80)
        for i in 1:80
            line[i] = read(io, Char)
        end
        push!(title, String(line))
    end
    skip(io, 4) # skip block size part

    # read header third block
    skip(io, 4) # skip block size part
    number_of_atom     = read(io, Int32)
    skip(io, 4) # skip block size part

    header_size = position(io)

    return DCDMetaInfo(file_size, first_step, nstep_save, total_step, total_unit,
                       time_step, unitcell_flag, version, number_of_atom, header_size)
end

"""
    read_dcd(filename::String;
             frame_indices::Union{Vector, OrdinalRange, Colon, Nothing} = nothing,
             step::Union{Int, Nothing} = nothing)
    ::Trajectory

Return Trajectory object which filled `coordinates`, `nframe`, `natom` fields.
If you set `frame_indices`, only specified frame is read.
If you set `step`, every step frame is read.
You can't specify both `frame_indices` and `step` at the same time.
"""
function read_dcd(filename::String;
                  frame_indices::Union{Vector, OrdinalRange, Colon, Nothing} = nothing,
                  step::Union{Int, Nothing} = nothing)::Trajectory
    coordinates_time_series = Matrix{Coordinate{Float32}}(undef, 0, 0)
    boxes_time_series = Vector{Coordinate{Float32}}()
    target_frame_indices = frame_indices

    open(filename, "r") do io

        dcd_meta::DCDMetaInfo = read_dcd_meta_info(io)
        header_size           = dcd_meta.header_size
        unitcell_flag         = dcd_meta.unitcell_flag
        number_of_atom        = dcd_meta.number_of_atom
        file_size             = dcd_meta.file_size

        # read body block
        seek(io, header_size)
        stepblocksize = (8 + 4 * number_of_atom) * 3
        if unitcell_flag
            stepblocksize += 56 # add unitcell block size (4 + 8 * 6 + 4)
        end
        total_frame = Int32(floor((file_size - header_size) / stepblocksize))
        if frame_indices != nothing && step != nothing
            throw(AssertionError("""
                                 You can't use both frame_indices and step option at the same time.
                                 """))
        elseif frame_indices != nothing
            if typeof(target_frame_indices) == Colon
                target_frame_indices = 1:total_frame
            else
                target_frame_indices = frame_indices
            end
        elseif step != nothing
            target_frame_indices = 1:step:total_frame
        else # both frame_indices and step are nothing
            target_frame_indices = 1:total_frame
        end

        coordinates_time_series =
            Matrix{Coordinate{Float32}}(undef, number_of_atom, length(target_frame_indices))
        output_frame_idx = 1
        for frame_idx in 1:total_frame
            # read unitcell block
            if unitcell_flag
                skip(io, 4) # skip block size part
                unitcell_info = Array{Float64, 1}(undef, 6)
                read!(io, unitcell_info)
                skip(io, 4) # skip block size part
                push!(boxes_time_series, Coordinate(unitcell_info[1],
                                                    unitcell_info[3],
                                                    unitcell_info[6]))
            end

            # read x coordinates
            skip(io, 4) # skip block size part
            x_coords = Array{Float32, 1}(undef, number_of_atom)
            read!(io, x_coords)
            skip(io, 4) # skip block size part

            # read y coorinates
            skip(io, 4) # skip block size part
            y_coords = Array{Float32, 1}(undef, number_of_atom)
            read!(io, y_coords)
            skip(io, 4) # skip block size part

            # read z coordinates
            skip(io, 4) # skip block size part
            z_coords = Array{Float32, 1}(undef ,number_of_atom)
            read!(io, z_coords)
            skip(io, 4) # skip block size part

            if frame_idx ∈ target_frame_indices
                atoms = collect(eachrow(hcat(x_coords, y_coords, z_coords)))
                coordinates_time_series[:, output_frame_idx] = Coordinate.(atoms)
                output_frame_idx += 1
            end
        end
    end
    Trajectory(coordinates_time_series, boxes=boxes_time_series)
end

"""
    write_dcd(filename::String, trj::Trajectory;
              save_step::Integer = 1, total_step::Integer = trj.nframe,
              unit_num::Integer  = 1,  time_step::Real = 1.0f0)

Write `coordinates`, `nframe` and `natom` information of `trj` to dcd file, named `filename`.
dcd file contain other information - save interval step, unit(chain) number in the system, total step and time step of original trajectory -, and if you do not specify these parameter, these are filled with default value.
"""
function write_dcd(filename::String, trj::Trajectory;
    save_step::Integer = 1, total_step::Integer = trj.nframe,
    unit_num::Integer = 1, time_step::Real = 1.0f0)

    open(filename, "w") do io

        # write header first block
        write(io, Int32(84))         # block size
        write(io, "CORD")            # signature
        write(io, Int32(trj.nframe)) # number of frames
        write(io, Int32(0))          # first frame index
        write(io, Int32(save_step))  # nstep_save
        write(io, Int32(total_step)) # number of steps
        write(io, Int32(unit_num))   # number of units
        write(io, zeros(Int32, 4))
        write(io, time_step)         # time step
        unitcell_flag = !isempty(trj.boxes)
        write(io, Int32(unitcell_flag))
        write(io, zeros(Int32, 8))
        write(io, Int32(24))         # version
        write(io, Int32(84))         # block size

        # write header second block
        write(io, Int32(84))         # block size
        write(io, Int32(1))          # number of lines
        write(io, "Generated by MolHandler.jl ((c) Yutaka Murata 2020) at $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))      ")
        write(io, Int32(84))

        # write header third block
        write(io, Int32(4))
        write(io, Int32(trj.natom))
        write(io, Int32(4))

        # write body block
        coord_part_size = 4 * trj.natom
        each_block_size = 8 + coord_part_size
        box_block_size  = 8 * 6
        for frame_idx in 1:trj.nframe

            # write box info block
            if unitcell_flag
                box_size = trj.boxes[frame_idx]
                write(io, Int32(box_block_size))
                write(io, Float64(box_size.x))
                write(io, Float64(90.0))
                write(io, Float64(box_size.y))
                write(io, Float64(90.0))
                write(io, Float64(90.0))
                write(io, Float64(box_size.z))
                write(io, Int32(box_block_size))
            end

            # write x coordinates
            write(io, Int32(coord_part_size))
            for atom_idx in 1:trj.natom
                write(io, Float32(trj.coordinates[atom_idx, frame_idx].x))
            end
            write(io, Int32(coord_part_size))

            # write y coordinates
            write(io, Int32(coord_part_size))
            for atom_idx in 1:trj.natom
                write(io, Float32(trj.coordinates[atom_idx, frame_idx].y))
            end
            write(io, Int32(coord_part_size))

            # write z coordinates
            write(io, Int32(coord_part_size))
            for atom_idx in 1:trj.natom
                write(io, Float32(trj.coordinates[atom_idx, frame_idx].z))
            end
            write(io, Int32(coord_part_size))
        end
    end
end

"""
    read_pdb(filename::AbstractString; model = :unspecified)::Trajectory

Return Trajectory object which filled all field.
If you set `:AA` to model field, this function set the mass of particle based on atomname field.
If you set `:CA` to model field, this function only read CA atom and set the mass of particcle mase on resname field.
"""
function read_pdb(filename::AbstractString; model = :unspecified)::Trajectory
    lines = open(filename, "r") do io
        readlines(io)
    end

    attributes = Vector{Attribute}()
    coordinates = Vector{Coordinate{Float32}}()
    first_frame = true
    if model == :unspecified || model == :AA
        for line in lines
            if occursin(r"^(ATOM|HETATM)", line)
                x_coord  = parse(Float32, line[31:38])
                y_coord  = parse(Float32, line[39:46])
                z_coord  = parse(Float32, line[47:54])
                push!(coordinates, Coordinate(x_coord, y_coord, z_coord))
                if first_frame
                    atomid   = parse(Int64, line[7:11])
                    atomname = strip(line[13:16])
                    resname  = strip(line[18:20])
                    resid    = parse(Int64, line[23:26])

                    if model == :unspecified
                        push!(attributes, Attribute(resname = resname, resid = resid,
                                                    atomname = atomname, atomid = atomid))
                    elseif model == :AA
                        push!(attributes, Attribute(resname = resname, resid = resid,
                                                    atomname = atomname, atomid = atomid,
                                                    mass = atom_mass(atomname)))
                    end
                end
            elseif occursin(r"^END", line) && first_frame
                first_frame = false
            end
        end
    elseif model == :CA
        for line in lines
            if occursin(r"^(ATOM|HETATM)", line) && strip(line[13:16]) == "CA"
                x_coord  = parse(Float32, line[31:38])
                y_coord  = parse(Float32, line[39:46])
                z_coord  = parse(Float32, line[47:54])
                push!(coordinates, Coordinate(x_coord, y_coord, z_coord))
                if first_frame
                    atomid   = parse(Int64, line[7:11])
                    resname  = strip(line[18:20])
                    resid    = parse(Int64, line[23:26])

                    push!(attributes, Attribute(resname = resname, resid = resid,
                                                atomname = "CA", atomid = atomid,
                                                mass = residue_mass(resname)))
                end
            elseif occursin(r"^END", line) && first_frame
                first_frame = false
            end
        end
    end

    atom_num  = length(attributes)
    @assert mod(length(coordinates), atom_num) == 0 "All frame should contain same number atoms."
    frame_num = div(length(coordinates), atom_num)
    Trajectory(reshape(coordinates, (atom_num, frame_num)), attributes=attributes)
end

"""
    write_pdb(filename::AbstractString, trj::Trajectory;
              tempfactor::Union{RealT, Nothing} = nothing) where RealT <: Real

Write `coordinates`, `resname`, `resid`, `atomname`, `atomid` to `filename` based on `trj`.
If you set Real value to `tempfactor`, all tempfactor will be field with that value.
"""
function write_pdb(filename::AbstractString, trj::Trajectory;
    tempfactor::Union{RealT, Nothing} = nothing,
    conects::Union{Vector{Vector{Int64}}, Nothing} = nothing) where RealT <: Real

    attributes = trj.attributes
    #  check attributes of trj have sufficient information.
    if any(attr->attr.resname==nothing, attributes)
        throw(AssertionError("""
                             There is nothing element in resname field of attribute in trajectory.
                             pdb format need resname for all particle.
                             """))
    elseif any(attr->attr.resid==nothing, attributes)
        throw(AssertionError("""
                             There is nothing element in resid field of attribute in trajectory.
                             pdb format need resid for all particle.
                             """))
    elseif any(attr->attr.atomname==nothing, attributes)
        throw(AssertionError("""
                             There is nothing element in atomname field of attribute in trajectory.
                             pdb format need atomname for all particle.
                             """))
    elseif any(attr->attr.atomid==nothing, attributes)
        throw(AssertionError("""
                             There is nothing element in atomid field of attribute in trajectory.
                             pdb format need atomid for all particle.
                             """))
    end

    coordinates = trj.coordinates
    # check attributes length match to the number of coordinates row.
    if (length(attributes) != size(coordinates)[1])
        throw(AssertionError("""
                             There is missmatch between the number of attribute and the nubmer of atom in trajectory.
                             """))
    end

    open(filename, "w") do io
        # write header
        write(io, "REMARK    Generated by MolHandler.jl ((c) Yutaka Murata 2020) at $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))\n")

        tempfactor_str = "      "
        if tempfactor != nothing
            tempfactor_str = Printf.@sprintf("%6.2f", tempfactor)
        end
        for (step, frame) in enumerate(eachcol(trj.coordinates))
            # write frame info
            write(io, "TITLE     step= $(step)\n")
            for (coordinate, attribute) in zip(frame, attributes)
                atomid_str   = Printf.@sprintf("%5d",   attribute.atomid)
                atomname_str = Printf.@sprintf("%-4s",  attribute.atomname)
                resname_str  = Printf.@sprintf("%3s",   attribute.resname)
                resid_str    = Printf.@sprintf("%4d",   attribute.resid)
                x_coord_str  = Printf.@sprintf("%8.3f", coordinate.x)
                y_coord_str  = Printf.@sprintf("%8.3f", coordinate.y)
                z_coord_str  = Printf.@sprintf("%8.3f", coordinate.z)
                write(io,
                      "ATOM  "*atomid_str*" "*atomname_str*" "*resname_str*"  "*resid_str*"    "*
                      x_coord_str*y_coord_str*z_coord_str*"      "*tempfactor_str*"\n")
            end
            write(io, "TER\n")

            # write CONECT records
            if conects != nothing
                for conect_info in conects
                    conect_str = "CONECT"
                    for idx in conect_info
                        conect_str *= Printf.@sprintf("%5d", idx)
                    end
                    conect_str *= "\n"
                    write(io, conect_str)
                end
            end
            write(io, "END\n")
        end
    end
end

"""
    read_xyz(filename::AbstractString)::Trajectory

Return Trajectory object which filled `coordinates`, `nframe`, `natom` fields.
"""
function read_xyz(filename::AbstractString)::Trajectory

    coordinates_time_series = Matrix{Coordinate{Float32}}(undef, 0, 0)
    attributes              = Vector{Attribute}()

    open(filename, "r") do fp
        # count number of atom
        number_of_atom = parse(Int, readline(fp))
        section_size = number_of_atom + 2

        # read atomname info
        readline(fp) # skip comment line
        for atom_idx in 1:number_of_atom
            atomname = split(readline(fp))[1]
            push!(attributes, Attribute(atomname = atomname))
        end

        # count number of frame
        seekstart(fp)
        line_num = countlines(fp)
        number_of_frame = div(line_num, section_size)

        # read coordinate
        coordinates_time_series =
            Matrix{Coordinate{Float32}}(undef, number_of_atom, number_of_frame)

        seekstart(fp)
        for frame_idx in 1:number_of_frame
            # skip number of atom and comment line
            for i in 1:2
                readline(fp)
            end
            for atom_idx in 1:number_of_atom
                line_elems = split(readline(fp))
                coord = Coordinate(map(c->parse(Float32, c), view(line_elems, 2:4)))
                coordinates_time_series[atom_idx, frame_idx] = coord
            end
        end
    end
    Trajectory(coordinates_time_series, attributes=attributes)
end

"""
    write_xyz(filename:AbstractString, trj::Trajectory)

Write `coordinates`, `nframe` and `natom` information of `trj` to xyz file, named `filename`.
"""
function write_xyz(filename::AbstractString, trj::Trajectory)
    atomname_arr = Vector{String}()
    for attr in trj.attributes
        if attr.atomname == nothing
            push!(atomname_arr, "UNK")
        else
            push!(atomname_arr, attr.atomname)
        end
    end

    natom = trj.natom
    open(filename, "w") do io
        for frame_idx in 1:trj.nframe
            write(io, string(natom)*"\n")
            write(io, "Generated by MolHandler.jl ((c) Yutaka Murata 2020)\n")
            for (atomname, coord) in zip(atomname_arr, view(trj.coordinates, :, frame_idx))
                x_coord_str = Printf.@sprintf("%.7f", coord.x)
                y_coord_str = Printf.@sprintf("%.7f", coord.y)
                z_coord_str = Printf.@sprintf("%.7f", coord.z)
                write(io, atomname*" "*x_coord_str*" "*y_coord_str*" "*z_coord_str*"\n")
            end
        end
    end
end
