function readdcd(filename::String)

    coordinates_time_series = Matrix{Atom}(undef, 0, 0)

    open(filename, "r") do io
        # read header first block
        skip(io, 4) # skip block size part
        header_sig = Array{Char, 1}(undef, 4)
        for i in 1:4
            header_sig[i] = read(io, Char)
        end
        total_frame = read(io, Int32)
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
#        println("header sig ", header_sig)
#        println("total frame ", total_frame)
#        println("first step ", first_step)
#        println("nstep save ", nstep_save)
#        println("total step ", total_step)
#        println("total unit ", total_unit)
#        println("time_step ", time_step)
#        println("version ", version)

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
#        println("title \n", title)

        # read header third block
        skip(io, 4) # skip block size part
        number_of_atom     = read(io, Int32)
        skip(io, 4) # skip block size part
#       println("number_of_atom ", number_of_atom)

        # read body block
        coordinates_time_series = Matrix{Atom}(undef, number_of_atom, total_frame)
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

            atoms = map(xyz -> Atom(xyz), eachrow(hcat(x_coords, y_coords, z_coords)))
            coordinates_time_series[:, frame_idx] = atoms
        end
    end
    Trajectory(coordinates_time_series)
end
