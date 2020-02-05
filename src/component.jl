struct Atom
    coordinate::Vector{Float64}
    resname::String
    resid::Int64
    atomname::String
    atomid::Int64
    mass::Float64
end

struct Frame
    atoms::Vector{Atom}
    natom::Int64
    nframe::Int64
end
