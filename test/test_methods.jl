@testset "readdcd function" begin
    trj = readdcd("data/test_position.dcd")
    @test isapprox(trj.coordinates[1,1],    [15.308, 14.180, -2.955], atol = 1e-3)
    @test isapprox(trj.coordinates[5,1001], [-2.577, 88.384, -7.513], atol = 1e-3)
end

@testset "readpdb function" begin
    trj = readpdb("data/1aki.pdb")
    @test isapprox(trj.coordinates[1,1], [35.365, 22.342, -11.980], atol = 1e-3)
    @test isapprox(trj.coordinates[1079, 1], [43.755, 23.843, 8.038], atol = 1e-3)
    @test trj.attributes[123].atomname == "N"
    @test trj.attributes[200].resid == 26
    @test trj.attributes[231].resname == "VAL"
end

@testset "get_frame" begin
    coordinates = Matrix{Vector{Float32}}(undef, 2, 2)
    coordinates[1,1] = [1.1f0, 1.2f0, 1.3f0]
    coordinates[1,2] = [2.1f0, 2.2f0, 2.3f0]
    coordinates[2,1] = [3.1f0, 3.2f0, 3.2f0]
    coordinates[2,2] = [4.1f0, 4.2f0, 4.1f0]
    attributes  = [Attribute(atomname = "hoge"), Attribute(atomname = "huga")]
    trj = Trajectory(coordinates, attributes)
    frame1 = get_frame(1, trj)
    frame2 = get_frame(2, trj)
    @test isapprox(frame1.coordinates[1][1], 1.1f0, atol = 1e-7)
    @test frame2.attributes[2].atomname   == "huga"
end

@testset "get_atom" begin
    coordinates = Matrix{Vector{Float32}}(undef, 2, 2)
    coordinates[1,1] = [1.1f0, 1.2f0, 1.3f0]
    coordinates[1,2] = [2.1f0, 2.2f0, 2.3f0]
    coordinates[2,1] = [3.1f0, 3.2f0, 3.2f0]
    coordinates[2,2] = [4.1f0, 4.2f0, 4.1f0]
    attributes  = [Attribute(atomname = "hoge"), Attribute(atomname = "huga")]
    trj = Trajectory(coordinates, attributes)
    atom1_series = get_atom(1, trj)
    atom2_series = get_atom(2, trj)
    @test isapprox(atom1_series[1].coordinate[1], 1.1f0, atol=1e-7)
    @test atom2_series[2].attribute.atomname == "huga"

    frame = Frame(coordinates[:,1], attributes)
    atom1 = get_atom(1, frame)
    atom2 = get_atom(2, frame)
    @test isapprox(atom1.coordinate[1], 1.1f0, atol = 1e-7)
    @test atom2.attribute.atomname == "huga"
end
