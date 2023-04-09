using LinearAlgebra,SpinFRGLattices
using LinearAlgebra:norm
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

abstract type AbstractSymop end

struct SiteTransformation2{T<:Real} <: AbstractSymop
    origin::SVector{3,T}
    matrix::SMatrix{3,3,T,9}
    WMatrix::SMatrix{4,4,T,16}
end

SiteTransformation2(origin::AbstractVector{T},matrix::AbstractMatrix{T}) where T <: Real = SiteTransformation2(origin,matrix,SiteTransformationMatrix(origin,matrix)) 

SiteTransformation = SiteTransformation2

function SiteTransformationMatrix(origin::AbstractVector{T},matrix::AbstractMatrix{T}) where T <: Real
    convert(SMatrix{4,4,T,16},[matrix origin; 0 0 0 1])
end

function (S::SiteTransformation2)(vec::AbstractVector{T}) where {T<:Real}
    x,y,z = vec
    x,y,z,_ = S.WMatrix * SA[x,y,z,one(x)]
    return SA[x,y,z]
end

parseAsSMatrix(v::AbstractVector{<:AbstractString}) = SMatrix{3,3,Float64,9}([parse(Float64,s) for s in v])'

function Base.show(io::IO, x::SiteTransformation2{T}) where T
    println(io,"3 dim SiteTransformation{T}")
    Base.print_matrix(io,x.WMatrix)
end
Base.show(io::IO, ::MIME"text/plain", x::AbstractSymop) = show(io,x)


"""read the symmetry operations from a vesta file"""
function readVestaSymops(filename)
    Info = readFileInfo(filename,"SYMOP",VestaFile())
    Info = getFields.(Info)
    origins = [parseAsSVector(data[1:3]) for data in Info]
    matrices = [parseAsSMatrix(data[4:12]) for data in Info]
    SiteTransformation.(origins,matrices)
end

"""for a position in fractional coordinates return true if it is in the unit cell"""
function isInUnitCell(pos_fractional::AbstractVector)
    all(pos_fractional .>= 0) && all(pos_fractional .< 1)
end

"""Given the symmetry inequivalent sites, the symmetry operations and the lattice vectors return the positions of all sites in a single unit cell by generating new sites via symmetry transformations until no new sites are found."""
function getUnitCell(sites::AbstractVector{<:AbstractVector},symops::AbstractVector{<:AbstractSymop},maxiter = 1000)
    positions = copy(sites)
    siteTypes = collect(eachindex(sites))
    i = 0
    for (x,site) in zip(siteTypes,positions)
        #new sites that are found are appended to positions and then used for generation iteratively. This loop should terminate eventually, since the number of sites in a unit cell is finite.
        #x is the type of the site, which does not change under symmetry transformations since, the site types are symmetry inequivalent by definition.
        for symop in symops
            r´ = symop(site)
            if isInUnitCell(r´) && r´ ∉ positions
                push!(positions, r´)
                push!(siteTypes, x)
            end
        end
        i +=1
        i >= maxiter && error("too many iterations")
    end
    perm = sortperm(positions)

    positions = positions[perm]
    siteTypes = siteTypes[perm]

    return (;positions,siteTypes)
end

function getUnitCell(filename)
    sites = readVestaSites(filename)
    symops = readVestaSymops(filename)
    getUnitCell(sites,symops)
end

function getBasis(filename)
    a1,a2,a3 = readLatticeVectors(filename)
    T = inv([a1 a2 a3])

    b,SiteType = getUnitCell(filename)
    
    b = [T * r for r in b]

    NUnique = maximum(SiteType)
    NNdist = min(getMinDistance(b),norm(a1),norm(a2),norm(a3))

    return Basis_Struct_3D(;a1,a2,a3,b,NUnique,SiteType,NNdist)
    # return (;a1,a2,a3,b,NUnique,SiteType,NNdist)
end

"""given a vector of sites return the minimum distance between any two sites"""
function getMinDistance(sites::AbstractVector{<:AbstractVector})
    minDist = Inf
    for i in eachindex(sites)
        for j in i+1:length(sites)
            minDist = min(minDist, norm(sites[i] - sites[j]))
        end
    end
    return minDist
end