import Dates
import Printf

"""
    read_dcd(filename::String;
             frame_indices::Union{Vector, OrdinalRange, Colon} = :)
    ::Trajectory

Return Trajectory object which filled `coordinates`, `nframe`, `natom` fields.
If you set frame_indices, only specified frame is read.
"""
function read_dcd(filename::String; frame_indices::Union{Vector, OrdinalRange, Colon} = :)::Trajectory
    coordinates_time_series = Matrix{Coordinate{Float32}}(undef, 0, 0)
    target_frame_indices = frame_indices

    open(filename, "r") do io
        unitcell_flag = false # for charmm format

        seekend(io)
        file_size = position(io) # get file size
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
        total_step     = read(io, Int32)
        total_unit     = read(io, Int32)
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

        # read body block
        header_size = position(io)
        stepblocksize = (8 + 4 * number_of_atom) * 3
        if unitcell_flag
            stepblocksize += 56 # add unitcell block size (4 + 8 * 6 + 4)
        end
        total_frame = Int32(floor((file_size - header_size) / stepblocksize))
        if typeof(target_frame_indices) == Colon
            target_frame_indices = 1:total_frame
        end

        coordinates_time_series =
            Matrix{Coordinate{Float32}}(undef, number_of_atom, length(target_frame_indices))
        output_frame_idx = 1
        for frame_idx in 1:total_frame
            # skip unitcell block
            if unitcell_flag
                skip(io, 56) # skip unitcell block size (4 + 8 * 6 + 4)
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
    Trajectory(coordinates_time_series)
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
        write(io, zeros(Int32, 9))
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
        for frame_idx in 1:trj.nframe
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
If you set `:CA` to model field, this function set the mass of particcle mase on resname field.
"""
function read_pdb(filename::AbstractString; model = :unspecified)::Trajectory
    lines = open(filename, "r") do io
        readlines(io)
    end

    attributes = Vector{Attribute}()
    coordinates = Vector{Coordinate{Float32}}()
    for line in lines
        if occursin(r"^(ATOM|HETATM)", line)
            atomid   = parse(Int64, line[7:11])
            atomname = strip(line[13:16])
            resname  = strip(line[18:20])
            resid    = parse(Int64, line[23:26])
            x_coord  = parse(Float32, line[31:38])
            y_coord  = parse(Float32, line[39:46])
            z_coord  = parse(Float32, line[47:54])

            if model == :unspecified
                push!(attributes, Attribute(resname = resname, resid = resid,
                                            atomname = atomname, atomid = atomid))
            elseif model == :AA
                push!(attributes, Attribute(resname = resname, resid = resid,
                                            atomname = atomname, atomid = atomid,
                                            mass = atom_mass(atomname)))
            elseif model == :CA
                push!(attributes, Attribute(resname = resname, resid = resid,
                                           atomname = atomname, atomid = atomid,
                                           mass = residue_mass(resname)))
            end
            push!(coordinates, Coordinate([x_coord, y_coord, z_coord]))
        end
    end
    Trajectory(reshape(coordinates, (length(coordinates), 1) ), attributes)
end

"""
    write_pdb(filename::AbstractString, trj::Trajectory;
              tempfactor::Union{RealT, Nothing} = nothing) where RealT <: Real

Write `coordinates`, `resname`, `resid`, `atomname`, `atomid` to `filename` based on `trj`.
If you set Real value to `tempfactor`, all tempfactor will be field with that value.
"""
function write_pdb(filename::AbstractString, trj::Trajectory;
                   tempfactor::Union{RealT, Nothing} = nothing) where RealT <: Real
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
            write(io, "ENDMDL\n")
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
    Trajectory(coordinates_time_series, attributes)
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
