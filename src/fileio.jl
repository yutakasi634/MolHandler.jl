function readdcd(filename::String)
    open(filename, "r") do io
        # read first block
        block_size_1 = read(io, Int32)
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
        block_size_2 = read(io, Int32)
        println("header sig ", header_sig)
        println("total frame ", total_frame)
        println("first step ", first_step)
        println("nstep save ", nstep_save)
        println("total step ", total_step)
        println("total unit ", total_unit)
        println("time_step ", time_step)
        println("version ", version)
        println("block size 2 ", block_size_2)
    end
end
