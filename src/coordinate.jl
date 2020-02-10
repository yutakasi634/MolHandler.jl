import Base: +, -, *, /, length, iterate
mutable struct Coordinate{T<:Real}
    # This class is for make vector operation broadcastable.
    x::T
    y::T
    z::T

    function Coordinate(vector::Array{T}) where T <: Real
        if length(vector) != 3
            ArgumentError("Argument for Coordinate constructor must be array of length 3")
        end
        new{T}(vector[1], vector[2], vector[3])
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

function *(arg1::Coordinate{T}, arg2::Coordinate{U}) where T <: Real where U <: Real
    arg1.x * arg2.x + arg1.y * arg2.y + arg1.z * arg2.z
end

function /(arg1::Coordinate{T}, arg2::U) where T <: Real where U <: Real
    Coordinate([arg1.x / arg2, arg1.y / arg2, arg1.z / arg2])
end

function length(arg::Coordinate{T}) where T <: Real
    1
end

function iterate(arg::Coordinate{T}) where T <: Real
    (arg, nothing)
end

function iterate(arg::Coordinate{T}, nothing) where T <: Real
    nothing
end

function norm(arg::Coordinate{T}) where T <: Real
    sqrt(abs2(arg.x) + abs2(arg.y) + abs2(arg.z))
end

function Array(arg::Coordinate{T}) where T <: Real
    [arg.x, arg.y, arg.z]
end
