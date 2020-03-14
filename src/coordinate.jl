import Base: +, -, *, /, length, iterate, zero, Array
"""
This correspond to coordinate of atom.
This contain 3 fields like below.

- x <: Real
- y <: Real
- z <: Real

Some operators and functions were overloaded to make Coordinate objects broadcastable.
"""
mutable struct Coordinate{RealT <: Real}
    # This class is for make vector operation broadcastable.
    x::RealT
    y::RealT
    z::RealT

    function Coordinate(vector::VectorT) where VectorT <: AbstractArray{RealT} where RealT <: Real
        if length(vector) != 3
            ArgumentError("Argument for Coordinate constructor must be array of length 3")
        end
        new{RealT}(vector[1], vector[2], vector[3])
    end
end

function +(arg1::Coordinate{T}, arg2::Coordinate{U}) where T <: Real where U <: Real
    Coordinate([arg1.x + arg2.x, arg1.y + arg2.y, arg1.z + arg2.z])
end

function -(arg1::Coordinate{T}, arg2::Coordinate{U}) where T <: Real where U <: Real
    Coordinate([arg1.x - arg2.x, arg1.y - arg2.y, arg1.z - arg2.z])
end

function *(arg1::Coordinate{T}, arg2::U) where T <: Real where U <: Real
    Coordinate([arg1.x * arg2, arg1.y * arg2, arg1.z * arg2])
end

function *(arg1::U, arg2::Coordinate{T}) where T <: Real where U <: Real
    Coordinate([arg1 * arg2.x, arg1 * arg2.y, arg1 * arg2.z])
end

function *(arg1::Coordinate{T}, arg2::Coordinate{U}) where T <: Real where U <: Real
    arg1.x * arg2.x + arg1.y * arg2.y + arg1.z * arg2.z
end

function /(arg1::Coordinate{T}, arg2::U) where T <: Real where U <: Real
    Coordinate([arg1.x / arg2, arg1.y / arg2, arg1.z / arg2])
end

# for operator broadcast
function length(arg::Coordinate{T}) where T <: Real
    1
end

function iterate(arg::Coordinate{T}) where T <: Real
    (arg, nothing)
end

function iterate(arg::Coordinate{T}, nothing) where T <: Real
    nothing
end

function zero(arg::Coordinate{T}) where T <: Real
    Coordinate([T(0.0) T(0.0) T(0.0)])
end

function Array(arg::Coordinate{T}) where T <: Real
    [arg.x, arg.y, arg.z]
end

function norm(arg::Coordinate{T}) where T <: Real
    sqrt(abs2(arg.x) + abs2(arg.y) + abs2(arg.z))
end

function distance(first_coord::Coordinate{T}, second_coord::Coordinate{T}) where T <: Real
    norm(first_coord - second_coord)
end
