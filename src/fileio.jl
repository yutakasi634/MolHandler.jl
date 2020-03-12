"""
    readdcd(filename::String)::Trajectory

Return Trajectory object which filled coordinates, `nframe`, `natom` fields.
"""
function readdcd(filename::String; frame_indices::Union{Vector, OrdinalRange, Colon} = :)::Trajectory

    coordinates_time_series = Matrix{Coordinate{Float32}}(undef, 0, 0)
    target_frame_indices = frame_indices

    open(filename, "r") do io
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
        coordblocksize = (8 + 4 * number_of_atom) * 3
        total_frame = Int32((file_size - header_size) / coordblocksize)
        if typeof(target_frame_indices) == Colon
            target_frame_indices = 1:total_frame
        end

        coordinates_time_series =
            Matrix{Coordinate{Float32}}(undef, number_of_atom, length(target_frame_indices))
        output_frame_idx = 1
        for frame_idx in 1:total_frame
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
    readpdb(filename::String)::Trajectory

Return Trajectory object which filled all field.
"""
function readpdb(filename::String)::Trajectory
    lines = open(filename, "r") do io
        readlines(io)
    end

    attributes = []
    coordinates = []
    for line in lines
        if occursin(r"^(ATOM|HETATM)", line)
            atomid   = parse(Int64, line[7:11])
            atomname = strip(line[13:16])
            resname  = strip(line[18:20])
            resid    = parse(Int64, line[23:26])
            x_coord  = parse(Float32, line[31:38])
            y_coord  = parse(Float32, line[39:46])
            z_coord  = parse(Float32, line[47:54])

            push!(attributes, Attribute(resname = resname, resid = resid,
                                        atomname = atomname, atomid = atomid))
            push!(coordinates, Coordinate([x_coord, y_coord, z_coord]))
        end
    end
    Trajectory(reshape(coordinates, (length(coordinates), 1) ), attributes)
end
