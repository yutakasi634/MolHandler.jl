function prepare_coordinates(frame_num::Int64)::Matrix{Coordinate{Float32}}
    coordinates = Matrix{Coordinate{Float32}}(undef, 3, frame_num)
    for atom_idx in 1:3, frame_idx in 1:frame_num
        coordinates[atom_idx, frame_idx] =
            Coordinate([atom_idx * 10.0f0 + frame_idx * 1.0f0 + 0.1f0,
                        atom_idx * 10.0f0 + frame_idx * 1.0f0 + 0.2f0,
                        atom_idx * 10.0f0 + frame_idx * 1.0f0 + 0.3f0])
    end
    coordinates
end

@testset "read_dcd function" begin
    trj = read_dcd("data/test_position.dcd")
    @test isapprox(Array(trj.coordinates[1,1]),    [15.308, 14.180, -2.955], atol = 1e-3)
    @test isapprox(Array(trj.coordinates[5,1001]), [-2.577, 88.384, -7.513], atol = 1e-3)

    trj_step = read_dcd("data/test_position.dcd", frame_indices = 1:10:1001)
    @test isapprox(Array(trj_step.coordinates[1,  1]), Array(trj.coordinates[1,    1]), atol = 1e-3)
    @test isapprox(Array(trj_step.coordinates[5,101]), Array(trj.coordinates[5, 1001]), atol = 1e-3)
end

@testset "write_dcd function" begin
    trj = read_dcd("data/test_position.dcd")
    write_dcd("data/write_test.dcd", trj)
    writed_trj = read_dcd("data/write_test.dcd")
    @test trj.nframe == writed_trj.nframe
    @test trj.natom  == writed_trj.natom
    @test isapprox(trj.coordinates[3, 10].x, writed_trj.coordinates[3, 10].x, atol = 1e-3)
end

@testset "read_pdb function" begin
    trj = read_pdb("data/1aki.pdb")
    @test isapprox(Array(trj.coordinates[1,1]), [35.365, 22.342, -11.980], atol = 1e-3)
    @test isapprox(Array(trj.coordinates[1079, 1]), [43.755, 23.843, 8.038], atol = 1e-3)
    @test trj.attributes[123].atomname == "N"
    @test trj.attributes[200].resid == 26
    @test trj.attributes[231].resname == "VAL"

    trj = read_pdb("data/1aki.pdb", model = :AA)
    @test isapprox(Array(trj.coordinates[1,1]), [35.365, 22.342, -11.980], atol = 1e-3)
    @test isapprox(Array(trj.coordinates[1079, 1]), [43.755, 23.843, 8.038], atol = 1e-3)
    @test trj.attributes[123].atomname == "N"
    @test trj.attributes[200].resid == 26
    @test trj.attributes[231].resname == "VAL"
    @test isapprox(trj.attributes[123].mass, 14.0069, atol = 1e-3)

    trj = read_pdb("data/1aki.pdb", model = :CA)
    @test isapprox(Array(trj.coordinates[1,1]), [35.365, 22.342, -11.980], atol = 1e-3)
    @test isapprox(Array(trj.coordinates[1079, 1]), [43.755, 23.843, 8.038], atol = 1e-3)
    @test trj.attributes[123].atomname == "N"
    @test trj.attributes[200].resid == 26
    @test trj.attributes[231].resname == "VAL"
    @test isapprox(trj.attributes[231].mass, 99.06841, atol = 1e-3)

    # test for multi frame pdb case
    trj = read_pdb("data/test_position.pdb")
    @test trj.nframe == 1001
    @test trj.natom  == 5
    @test isapprox(Array(trj.coordinates[1,   1]), [15.308, 14.180, -2.955], atol = 1e-3)
    @test isapprox(Array(trj.coordinates[5,1001]), [-2.578, 88.385, -7.513] ,atol = 1e-3)
end

@testset "write_pdb function" begin
    trj = read_pdb("data/1aki.pdb")
    write_pdb("data/write_test.pdb", trj)
    writed_trj = read_pdb("data/write_test.pdb")
    @test isapprox(trj.coordinates[1, 1].x, writed_trj.coordinates[1, 1].x, atol = 1e-3)
    @test isapprox(trj.coordinates[1079, 1].x, writed_trj.coordinates[1079, 1].x, atol = 1e-3)
    @test trj.attributes[123].atomname == writed_trj.attributes[123].atomname
    @test trj.attributes[200].resid    == writed_trj.attributes[200].resid
    @test trj.attributes[231].resname  == writed_trj.attributes[231].resname

    write_pdb("data/write_test.pdb", trj, tempfactor = 0.00)
    @test isapprox(trj.coordinates[1, 1].x, writed_trj.coordinates[1, 1].x, atol = 1e-3)
    @test isapprox(trj.coordinates[1079, 1].x, writed_trj.coordinates[1079, 1].x, atol = 1e-3)
    @test trj.attributes[123].atomname == writed_trj.attributes[123].atomname
    @test trj.attributes[200].resid    == writed_trj.attributes[200].resid
    @test trj.attributes[231].resname  == writed_trj.attributes[231].resname

    temp_arr = []
    open("data/write_test.pdb", "r") do fp
        for line in eachline(fp)
            if occursin(r"^ATOM", line)
                push!(temp_arr, line[63:66])
            end
        end
    end
    @test temp_arr[5]   == "0.00"
    @test temp_arr[756] == "0.00"
end

@testset "read_xyz function" begin
    trj = read_xyz("data/test_position.xyz")
    @test isapprox(Array(trj.coordinates[1, 1]),    [15.308, 14.180, -2.955], atol = 1e-3)
    @test isapprox(Array(trj.coordinates[5, 1001]), [-2.577, 88.384, -7.513], atol = 1e-3)
    @test trj.attributes[3].atomname == "CA"
end

@testset "write_xyz function" begin
    trj = read_xyz("data/test_position.xyz")
    write_xyz("data/write_test.xyz", trj)
    writed_trj = read_xyz("data/write_test.xyz")
    @test trj.nframe == writed_trj.nframe
    @test trj.natom  == writed_trj.natom
    @test writed_trj.attributes[3].atomname == trj.attributes[3].atomname
    @test isapprox(trj.coordinates[3, 10].x, writed_trj.coordinates[3, 10].x, atol = 1e-3)

    trj = read_dcd("data/test_position.dcd")
    write_xyz("data/write_test.xyz", trj)
    writed_trj = read_xyz("data/write_test.xyz")
    @test writed_trj.attributes[3].atomname == "UNK"
end

@testset "get_frame" begin
    coordinates = Matrix{Coordinate{Float32}}(undef, 2, 2)
    coordinates[1,1] = Coordinate([1.1f0, 1.2f0, 1.3f0])
    coordinates[1,2] = Coordinate([2.1f0, 2.2f0, 2.3f0])
    coordinates[2,1] = Coordinate([3.1f0, 3.2f0, 3.2f0])
    coordinates[2,2] = Coordinate([4.1f0, 4.2f0, 4.1f0])
    attributes  = [Attribute(atomname = "hoge"), Attribute(atomname = "huga")]
    trj = Trajectory(coordinates, attributes)
    frame1 = get_frame(1, trj)
    frame2 = get_frame(2, trj)
    @test isapprox(Array(frame1.coordinates[1])[1], 1.1f0, atol = 1e-7)
    @test frame2.attributes[2].atomname   == "huga"
end

@testset "get_atom" begin
    coordinates = Matrix{Coordinate{Float32}}(undef, 2, 2)
    coordinates[1,1] = Coordinate([1.1f0, 1.2f0, 1.3f0])
    coordinates[1,2] = Coordinate([2.1f0, 2.2f0, 2.3f0])
    coordinates[2,1] = Coordinate([3.1f0, 3.2f0, 3.2f0])
    coordinates[2,2] = Coordinate([4.1f0, 4.2f0, 4.1f0])
    attributes  = [Attribute(atomname = "hoge"), Attribute(atomname = "huga")]
    trj = Trajectory(coordinates, attributes)
    atom1_series = get_atom(1, trj)
    atom2_series = get_atom(2, trj)
    @test isapprox(atom1_series[1].coordinate.x, 1.1f0, atol=1e-7)
    @test atom2_series[2].attribute.atomname == "huga"

    frame = Frame(coordinates[:,1], attributes)
    atom1 = get_atom(1, frame)
    atom2 = get_atom(2, frame)
    @test isapprox(atom1.coordinate.x, 1.1f0, atol = 1e-7)
    @test atom2.attribute.atomname == "huga"
end

@testset "clip_trajectory(query::Integer)" begin
    coordinates = prepare_coordinates(3)
    attributes = [Attribute(atomname = "hoge"), Attribute(atomname = "huga"), Attribute(atomname = "piyo")]
    trj = Trajectory(coordinates, attributes)

    cliped_frame_int = clip_trajectory(2, trj)
    @test isapprox(cliped_frame_int.coordinates[2, 1].x, trj.coordinates[2, 2].x, atol = 1e-3)
    @test cliped_frame_int.attributes   == trj.attributes

    cliped_atom_int = clip_trajectory(2, trj, query_key = :atom)
    @test isapprox(cliped_atom_int.coordinates[1, 2].x,  trj.coordinates[2, 2].x, atol = 1e-3)
    @test cliped_atom_int.attributes   == [trj.attributes[2]]

end

@testset "clip_trajectory(query::Union{Array, StepRange})" begin
    coordinates = prepare_coordinates(3)
    attributes = [Attribute(atomname = "hoge"), Attribute(atomname = "huga"), Attribute(atomname = "piyo")]
    trj = Trajectory(coordinates, attributes)

    cliped_frame_array = clip_trajectory([2, 3], trj)
    @test isapprox(cliped_frame_array.coordinates[2, 2].x, trj.coordinates[2, 3].x, atol = 1e-3)
    @test cliped_frame_array.attributes  == trj.attributes

    cliped_frame_range = clip_trajectory(1:2:3, trj)
    @test isapprox(cliped_frame_range.coordinates[2, 2].y, trj.coordinates[2, 3].y, atol = 1e-3)
    @test cliped_frame_range.attributes  == trj.attributes

    cliped_atom_array = clip_trajectory([2, 3], trj, query_key = :atom)
    @test isapprox(cliped_atom_array.coordinates[2, 2].z,  trj.coordinates[3, 2].z, atol = 1e-3)
    @test cliped_atom_array.attributes  == trj.attributes[[2, 3]]

    cliped_atom_range = clip_trajectory(1:2:3, trj, query_key = :atom)
    @test isapprox(cliped_atom_range.coordinates[2, 2].x,  trj.coordinates[3, 2].x, atol = 1e-3)
    @test cliped_atom_range.attributes  == trj.attributes[1:2:3]
end

@testset "center_of_mass" begin
    coordinates = prepare_coordinates(3)
    trj_wo_attr = Trajectory(coordinates)

    @test_throws ArgumentError center_of_mass(trj_wo_attr)

    attributes = [Attribute(mass = 2.0f0), Attribute(mass = 3.0f0), Attribute(mass = 4.0f0)]
    trj = Trajectory(coordinates, attributes)

    @test isapprox(Array(center_of_mass(trj)[1]), [23.3222f0, 23.4222f0, 23.5222f0], atol = 1e-3)
    @test isapprox(Array(center_of_mass(trj, atom_indices = 1:2:3)[2]), [25.4333f0, 25.5333f0, 25.6333f0],
                   atol = 1e-3)
    @test isapprox(Array(center_of_mass(trj, atom_indices = 1:2:3, frame_indices = 2:3)[1]),
                   [25.4333f0, 25.5333f0, 25.6333f0], atol = 1e-3)

    @test isapprox(Array(center_of_mass(trj, geometric = true)[1]), [21.1f0, 21.2f0, 21.3f0], atol = 1e-3)
    @test isapprox(Array(center_of_mass(trj, geometric = true, atom_indices = 1:2:3)[2]),
                   [22.1f0, 22.2f0, 22.3f0], atol = 1e-3)
    @test isapprox(Array(center_of_mass(trj, geometric = true, atom_indices = 1:2:3, frame_indices = 2:3)[1]),
                   [22.1f0, 22.2f0, 22.3f0], atol = 1e-3)

end

@testset "pair_length_matrix" begin
    coordinates = prepare_coordinates(3)
    length_matrix = pair_length_matrix(coordinates[:, 1], coordinates[:, 3])

    @test size(length_matrix) == (3, 3)
    @test isapprox(38.1051f0, length_matrix[1, 3], atol = 1e-3)

    trj = Trajectory(coordinates)
    length_matrices = pair_length_matrix(trj, frame_indices = 1:2,
                                         first_atom_indices = 1:2:3,
                                         second_atom_indices = :)

    @test size(length_matrices[1]) == (2, 3)
    @test isapprox(34.6410f0, length_matrices[2][2, 1], atol = 1e-3)
end

@testset "pair_length_matrix_pbc" begin
    coordinates = prepare_coordinates(3)
    upper_bound = Coordinate(40.0, 40.0, 40.0)
    lower_bound = Coordinate( 0.0,  0.0,  0.0)
    length_matrix = pair_length_matrix_pbc(coordinates[:, 1], coordinates[:, 3],
                                           upper_bound, lower_bound)

    @test size(length_matrix) == (3, 3)
    @test isapprox(31.1769, length_matrix[1, 3], atol = 1e-3)

    trj = Trajectory(coordinates)
    length_matrices = pair_length_matrix_pbc(trj, upper_bound, lower_bound,
                                             frame_indices = 1:2:3,
                                             first_atom_indices = 1:2:3,
                                             second_atom_indices = :)

    @test size(length_matrices[1]) == (2, 3)
    @test isapprox(distance_pbc(coordinates[3, 3], coordinates[1, 3], upper_bound, lower_bound),
                   length_matrices[2][2, 1], atol = 1e-3)
end

@testset "pair_length_matrix_parallel" begin
    coordinates = prepare_coordinates(3)
    trj = Trajectory(coordinates)

    @test isapprox(pair_length_matrix(trj), pair_length_matrix_parallel(trj), atol = 1e-3)
end

@testset "contact_bool_matrix" begin
    coordinates = prepare_coordinates(2)
    bool_matrix = contact_bool_matrix(20.0f0, coordinates[:, 1], coordinates[:, 2])

    @test size(bool_matrix) == (3, 3)
    @test bool_matrix == [true true false; true true true; false true true]

    trj = Trajectory(coordinates)
    contact_bool_mat_arr = contact_bool_matrix(18.0f0, trj, first_atom_indices = 1:2, second_atom_indices = 2:3)

    @test contact_bool_mat_arr[1][1, 1] == true
    @test contact_bool_mat_arr[2][1, 2] == false
end

@testset "contact_bool_matrix_pbc" begin
    first_coord  = Coordinate(1.0, 2.0, 2.0)
    second_coord = Coordinate(2.5, 2.0, 2.0)
    third_coord  = Coordinate(5.0, 2.0, 2.0)
    coordinates  = [first_coord  first_coord ;
                    second_coord second_coord;
                    third_coord  third_coord ]
    trj          = Trajectory(coordinates)
    upper_bound  = Coordinate(5.5, 4.0, 4.0)
    lower_bound  = Coordinate(0.0, 0.0, 0.0)

    trj = Trajectory(coordinates)
    contact_bool_mat_arr = contact_bool_matrix_pbc(2.0f0, trj, upper_bound, lower_bound)

    @test contact_bool_mat_arr[1][1, 2] == true
    @test contact_bool_mat_arr[2][3, 1] == true
    @test contact_bool_mat_arr[2][2, 3] == false
end


@testset "contact_bool_matrix_parallel" begin
    coordinates = prepare_coordinates(3)
    trj = Trajectory(coordinates)

    @test isapprox(contact_bool_matrix(18.0f0, trj), contact_bool_matrix_parallel(18.0f0, trj), atol = 1e-3)
end

@testset "contact_probability_matrix" begin
    coordinates = prepare_coordinates(3)
    trj = Trajectory(coordinates)
    contact_prob_mat = contact_probability_matrix(18.0f0, trj,
                                                  frame_indices = 1:2:3,
                                                  first_atom_indices = 1:2,
                                                  second_atom_indices = 2:3)

    @test size(contact_prob_mat) == (2, 2)
    @test isapprox(contact_prob_mat[1, 1], 1.0f0, atol = 1e-3)
    @test isapprox(contact_prob_mat[1, 2], 0.0f0, atol = 1e-3)
end

@testset "contact_probability_matrix_parallel" begin
    coordinates = prepare_coordinates(3)
    trj = Trajectory(coordinates)

    @test isapprox(contact_probability_matrix(18.0f0, trj), contact_probability_matrix_parallel(18.0f0, trj), atol = 1e-3)
end

@testset "radius_of_gyration" begin
    coordinates = prepare_coordinates(3)
    attributes = [Attribute(mass = 2.0f0), Attribute(mass = 3.0f0), Attribute(mass = 4.0f0)]
    trj = Trajectory(coordinates, attributes)

    rg_arr     = Array{Float32, 1}(undef, 3)
    geo_rg_arr = Array{Float32, 1}(undef, 3)
    for frame_idx in 1:size(coordinates, 2)
        # com calculation
        com     = Coordinate([0.0, 0.0, 0.0])
        geo_com = Coordinate([0.0, 0.0, 0.0])
        for atom_idx in 1:size(coordinates, 1)
            com     += coordinates[atom_idx, frame_idx] * attributes[atom_idx].mass
            geo_com += coordinates[atom_idx, frame_idx]
        end
        com     /= sum([attr.mass for attr in attributes])
        geo_com /= size(coordinates, 1)

        # rg calculation
        rg     = 0.0
        geo_rg = 0.0
        for atom_idx in 1:size(coordinates, 1)
            rg     += norm(coordinates[atom_idx, frame_idx] - com)^2 * attributes[atom_idx].mass
            geo_rg += norm(coordinates[atom_idx, frame_idx] - geo_com)^2
        end
        rg     = sqrt(rg     / sum([attr.mass for attr in attributes]))
        geo_rg = sqrt(geo_rg / size(coordinates, 1))
        rg_arr[frame_idx]     = rg
        geo_rg_arr[frame_idx] = geo_rg
    end

    @test isapprox(radius_of_gyration(trj, frame_indices = 1:2:3),                   rg_arr[1:2:3],     atol = 1e-3)
    @test isapprox(radius_of_gyration(trj, frame_indices = 1:2:3, geometric = true), geo_rg_arr[1:2:3], atol = 1e-3)
end

@testset "atom_mass" begin
    @test isapprox(atom_mass("HB3"), 1.00798, atol = 1e-5)
    @test isapprox(atom_mass("CA"),  12.0106, atol = 1e-5)
    @test isapprox(atom_mass("NB"),  14.0069, atol = 1e-5)
    @test isapprox(atom_mass("NA"),  22.9898, atol = 1e-5)
    @test isapprox(atom_mass("CL"),  35.4520, atol = 1e-5)
end

@testset "residue_mass" begin
    @test isapprox(residue_mass("ALA"), 71.03711, atol = 1e-5)
    @test isapprox(residue_mass("HOH"), 18.01528, atol = 1e-5)
end

@testset "distance_pbc" begin
    first_coord  = Coordinate(1.0, 2.0, 2.0)
    second_coord = Coordinate(2.5, 2.0, 2.0)
    third_coord  = Coordinate(3.5, 2.0, 2.0)
    upper_bound  = Coordinate(4.0, 4.0, 4.0)
    lower_bound  = Coordinate(0.0, 0.0, 0.0)
    first_second =
        distance_pbc(first_coord, second_coord, upper_bound, lower_bound)
    first_third  =
        distance_pbc(first_coord, third_coord,  upper_bound, lower_bound)
    @test isapprox(first_second, 1.5, atol = 1e-3)
    @test isapprox(first_third,  1.5, atol = 1e-3)
end

@testset "fix_pbc!" begin
    dcd = read_dcd("data/test_position.dcd")
    @test_nowarn fix_pbc!(dcd, Coordinate([0.0f0, 0.0f0, 0.0f0]), Coordinate([10.0f0, 10.0f0, 10.0f0]))
end

@testset "move_pbc_center" begin
    first_coord  = Coordinate(1.0, 2.0, 2.0)
    second_coord = Coordinate(2.5, 2.0, 2.0)
    third_coord  = Coordinate(3.5, 2.0, 2.0)
    coordinates  = [first_coord, second_coord, third_coord]
    new_center   = Coordinate(1.0, 2.0, 2.0)
    box_size     = Coordinate(4.0, 4.0, 4.0)
    new_coords = move_pbc_center(coordinates, new_center, box_size)
    @test isapprox(Array(first_coord), Array(new_coords[1]), atol = 1e-3)
    @test isapprox([-0.5, 2.0, 2.0],   Array(new_coords[3]), atol = 1e-3)
end
