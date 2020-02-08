using MolHandler
using Test

@testset "readdcd function" begin
    trj = readdcd("test_position.dcd")
    @test isapprox(trj.coordinates[1,1],    [15.308, 14.180, -2.955], atol = 1e-3)
    @test isapprox(trj.coordinates[5,1001], [-2.577, 88.384, -7.513], atol = 1e-3)
end
