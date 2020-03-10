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

@testset "readdcd function" begin
    trj = readdcd("data/test_position.dcd")
    @test isapprox(Array(trj.coordinates[1,1]),    [15.308, 14.180, -2.955], atol = 1e-3)
    @test isapprox(Array(trj.coordinates[5,1001]), [-2.577, 88.384, -7.513], atol = 1e-3)
end

@testset "readpdb function" begin
    trj = readpdb("data/1aki.pdb")
    @test isapprox(Array(trj.coordinates[1,1]), [35.365, 22.342, -11.980], atol = 1e-3)
    @test isapprox(Array(trj.coordinates[1079, 1]), [43.755, 23.843, 8.038], atol = 1e-3)
    @test trj.attributes[123].atomname == "N"
    @test trj.attributes[200].resid == 26
    @test trj.attributes[231].resname == "VAL"
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
    @test cliped_frame_int.coordinates  == reshape(trj.coordinates[:, 2], (3, 1))
    @test cliped_frame_int.attributes   == trj.attributes

    cliped_atom_int = clip_trajectory(2, trj, query_key = :atom)
    @test cliped_atom_int.coordinates  == reshape(trj.coordinates[2, :], (1, 3))
    @test cliped_atom_int.attributes   == [trj.attributes[2]]

end

@testset "clip_trajectory(query::Union{Array, StepRange})" begin
    coordinates = prepare_coordinates(3)
    attributes = [Attribute(atomname = "hoge"), Attribute(atomname = "huga"), Attribute(atomname = "piyo")]
    trj = Trajectory(coordinates, attributes)

    cliped_frame_array = clip_trajectory([2, 3], trj)
    @test cliped_frame_array.coordinates == trj.coordinates[:, [2, 3]]
    @test cliped_frame_array.attributes  == trj.attributes

    cliped_frame_range = clip_trajectory(1:2:3, trj)
    @test cliped_frame_range.coordinates == trj.coordinates[:, 1:2:3]
    @test cliped_frame_range.attributes  == trj.attributes

    cliped_atom_array = clip_trajectory([2, 3], trj, query_key = :atom)
    @test cliped_atom_array.coordinates == trj.coordinates[[2, 3], :]
    @test cliped_atom_array.attributes  == trj.attributes[[2, 3]]

    cliped_atom_range = clip_trajectory(1:2:3, trj, query_key = :atom)
    @test cliped_atom_range.coordinates == trj.coordinates[1:2:3, :]
    @test cliped_atom_range.attributes  == trj.attributes[1:2:3]
end

@testset "center_of_mass" begin
    coordinates = prepare_coordinates(3)
    attributes = [Attribute(mass = 2.0f0), Attribute(mass = 3.0f0), Attribute(mass = 4.0f0)]
    trj = Trajectory(coordinates, attributes)
    @test isapprox(Array(center_of_mass(trj)[1]), [23.3222f0, 23.4222f0, 23.5222f0], atol = 1e-3)
    @test isapprox(Array(center_of_mass(trj, indices = 1:2:3)[2]), [25.4333f0, 25.5333f0, 25.6333f0],
                   atol = 1e-3)

    @test isapprox(Array(center_of_mass(trj, geometric = true)[1]), [21.1f0, 21.2f0, 21.3f0], atol = 1e-3)
    @test isapprox(Array(center_of_mass(trj, geometric = true, indices = 1:2:3)[2]),
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
