using LinearAlgebra

function containsData(line::String,::VestaFile)
    line = replace(line," " => "")
    return !isempty(line) && isdigit(first(line))
end

"""read a vesta file and return the angles and lengths of the lattice vectors under CELLP"""
function readCellp(filename)
    Info = readFileInfo(filename,"CELLP",VestaFile())[1]
    Info = getFields.(Info)
    a,b,c, alpha,beta,gamma = parse.(Float64,Info)
    (;a,b,c,alpha,beta,gamma)
end

function readLatticeVectors(filename)
    a,b,c, alpha,beta,gamma = readCellp(filename)
    getLatticeVectors(a,b,c,alpha,beta,gamma)
end

"""given three angles α,β,γ return a set of unit vectors a = [1,0,0], b = [b1,b2,0], c = [c1,c2,c3], such that the α is the angle between b and c, β is the angle between a and c and γ is the angle between a and b"""
function getUnitVectors(α, β, γ)
    @assert γ !== 180 "γ = 180° is not allowed"
    a = SA[1,0,0]
    b = SA[cosd(γ), sind(γ), 0]
    c1 = cosd(β)
    c2 = (cosd(α) - cosd(β) * cosd(γ)) / sind(γ)
    c3 = sqrt(abs(1 - c1^2 - c2^2))
    c = SA[c1,c2,c3] |> normalize
    return (a, b, c)
end

function latticeparameters(a::AbstractVector,b::AbstractVector,c::AbstractVector)
    a,b,c = map(normalize,(a,b,c))
    
    (;α = acosd(b' * c), β = acosd(a' * c), γ = acosd(a' * b))
end

"""given the lengths a,b,c and the angles alpha,beta,gamma return the lattice vectors"""
function getLatticeVectors(a, b, c, alpha, beta, gamma)
    avec,bvec,cvec = getUnitVectors(alpha,beta,gamma)
    avec *= a
    bvec *= b
    cvec *= c
    return avec,bvec,cvec
end

"""reads position of symmetry inequivalent sites in fractional coordinates"""
function readVestaSites(filename)
    Info = readFileInfo(filename,"STRUC",VestaFile())[1:2:end-1] # skip standard uncertainty lines and line of trailing zeros
    Info = getFields.(Info)
    pos = [parseAsSVector(data[5:7]) for data in Info]
end
