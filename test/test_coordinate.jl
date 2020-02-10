@testset "Coordinate constructor" begin
    coord = Coordinate([1.0, 2.0, 3.0])
    @test isapprox(coord.x, 1.0, atol = 1e-5)
    @test isapprox(coord.y, 2.0, atol = 1e-5)
    @test isapprox(coord.z, 3.0, atol = 1e-5)
end

@testset "operator + for Coordinate" begin
    coord1 = Coordinate([1.0, 2.0, 3.0])
    coord2 = Coordinate([2.0, 3.0, 4.0])
    res = coord1 + coord2
    @test isapprox(res.x, 3.0, atol = 1e-5)
    @test isapprox(res.z, 7.0, atol = 1e-5)
end

@testset "operator - for Coordinate" begin
    coord1 = Coordinate([3.0, 2.0, 1.0])
    coord2 = Coordinate([2.0, 3.0, 4.0])
    res = coord2 - coord1
    @test isapprox(res.x, -1.0, atol = 1e-5)
    @test isapprox(res.z, 3.0, atol = 1e-5)
end

@testset "operator * for Coordinate" begin
    coord1 = Coordinate([1.0, 2.0, 3.0])
    coord2 = Coordinate([2.0, 3.0, 4.0])
    coord_coord_case = coord1 * coord2
    @test isapprox(coord_coord_case, 20.0, atol = 1e-5)
    coord_scholar_case = coord1 * 2.0
    @test isapprox(Array(coord_scholar_case), [2.0, 4.0, 6.0], atol = 1e-5)

end

@testset "operator / for Coordinate" begin
    coord1 = Coordinate([1.0, 2.0, 3.0])
    res = coord1 / 2.0
    @test isapprox(res.x, 0.5, atol = 1e-5)
    @test isapprox(res.z, 1.5, atol = 1e-5)
end

@testset "norm for Coordinate" begin
    coord = Coordinate([1.0, 2.0, 3.0])
    resnorm = norm(coord)
    predicted_val    = sqrt(1.0^2 + 2.0^2 + 3.0^2)
    @test isapprox(resnorm, predicted_val, atol = 1e-5)
end

@testset "broadcast for Coordinate" begin
    coord_a = Coordinate([1.0, 2.0, 3.0])
    coord_b = Coordinate([2.0, 3.0, 4.0])
    coord_c = Coordinate([3.0, 4.0, 5.0])
    coord_d = Coordinate([3.0, 2.0, 1.0])
    coords_vec = [coord_a, coord_b, coord_c]
    @test isapprox(Array((coords_vec .+ coord_d)[3]), [6.0, 6.0, 6.0],    atol = 1e-5)
    @test isapprox(Array((coords_vec .- coord_d)[2]), [-1.0, 1.0, 3.0],   atol = 1e-5)
    @test isapprox(       coords_vec .* coord_d,      [10.0, 16.0, 22.0], atol = 1e-5)
    @test isapprox(Array((coords_vec .* 2.0)[2]),     [4.0, 6.0, 8.0],    atol = 1e-5)
end
