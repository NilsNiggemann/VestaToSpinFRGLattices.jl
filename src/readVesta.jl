using LinearAlgebra
"""read a vesta file and return the angles and lengths of the lattice vectors under CELLP"""
function readCellp(filename)
    Info = readCifInfo(filename,"CELLP")[1]
    Info = getFields.(Info)
    a,b,c, alpha,beta,gamma = parse.(Float64,Info)
    (;a,b,c,alpha,beta,gamma)
end

"""given three angles α,β,γ return a set of unit vectors a = [1,0,0], b = [b1,b2,0], c = [c1,c2,c3], such that the α is the angle between b and c, β is the angle between a and c and γ is the angle between a and b"""
function getUnitVectors(α, β, γ)
    a = SA[1,0,0]
    b = SA[cos(γ), sin(γ), 0]
    c1 = cos(β)
    c2 = (cos(α) - cos(β) * cos(γ)) / sin(γ)
    c3 = sqrt(abs(1 - c1^2 - c2^2))
    c = SA[c1,c2,c3] |> normalize
    return (a, b, c)
end

function latticeparameters(a::AbstractVector,b::AbstractVector,c::AbstractVector)
    a,b,c = map(normalize,(a,b,c))
    
    (;α = acos(b' * c), β = acos(a' * c), γ = acos(a' * b))
end

"""given the lengths a,b,c and the angles alpha,beta,gamma return the lattice vectors"""
function getLatticeVectors(a, b, c, alpha, beta, gamma)
    avec,bvec,cvec = getUnitVectors(alpha,beta,gamma)
    avec *= a
    bvec *= b
    cvec *= c
    return avec,bvec,cvec
end