using LinearAlgebra,SpinFRGLattices
using LinearAlgebra:norm
function containsData(line::String,::VestaFile)
    line = replace(line," " => "")
    return !isempty(line) && isdigit(first(line)) || first(line) ∈ ('-','+')
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
    function parseSites(s)
        number,atom,label,occ,x,y,z,_ = getFields(s)
        position = parse.(Float64,SA[x,y,z])
        return (;atom,label,position)
    end
    return parseSites.(Info)
end


abstract type AbstractSymop <: Function end

struct SiteTransformation{T<:Real} <: AbstractSymop
    WMatrix::SMatrix{4,4,T,16}
end

getOrigin(S::SiteTransformation) = S.WMatrix[1:3,4]
getRotation(S::SiteTransformation) = S.WMatrix[1:3,1:3]

SiteTransformation(origin::AbstractVector{T},matrix::AbstractMatrix{T}) where T <: Real = SiteTransformation(convert(SMatrix{4,4,T,16},[matrix origin; 0 0 0 1])) 

function (S::SiteTransformation)(vec::AbstractVector{T}) where {T<:Real}
    x,y,z = vec
    x,y,z,_ = S.WMatrix * SA[x,y,z,one(x)]
    return SA[x,y,z]
end

parseAsSMatrix(v::AbstractVector{<:AbstractString}) = SMatrix{3,3,Float64,9}([parse(Float64,s) for s in v])'

function Base.show(io::IO, x::SiteTransformation{T}) where T
    println(io,"3 dim SiteTransformation{",T,"}")
    Base.print_matrix(io,x.WMatrix)
end
Base.show(io::IO, ::MIME"text/plain", x::SiteTransformation) = show(io,x)

ispuretranslation(S::SiteTransformation) = abs(det(getRotation(S))) == 0

"""read the symmetry operations from a vesta file"""
function readVestaSymops(filename)
    Info = readFileInfo(filename,"SYMOP",VestaFile())
    Info = getFields.(Info)
    origins = [parseAsSVector(data[1:3]) for data in Info]
    matrices = [parseAsSMatrix(data[4:12]) for data in Info]
    Syms = SiteTransformation.(origins,matrices)
    filter!(!ispuretranslation,Syms)
    return Syms
end

"""for a position in fractional coordinates return true if it is in the unit cell"""
function isInUnitCell(pos_fractional::AbstractVector)
    all(pos_fractional .>= -0) && all(pos_fractional .< 1)
end

translate_back_to_UC(r::AbstractVector) = r .- floor.(r)
translate_back_to_UC(R::Rvec_3D) = Rvec(0,0,0,R.b)

"""given a vector and a vector of vectors return if the vector is approximately equal to one of the vectors in the vector of vectors"""
function approxin(r::AbstractVector,rs::AbstractVector{<:AbstractVector})
    any(norm(r - rx) < 1e-9 for rx in rs)
end

"""Given the symmetry inequivalent sites, the symmetry operations and the lattice vectors return the positions of all sites in a single unit cell by generating new sites via symmetry transformations until no new sites are found."""
function getUnitCell(sites::AbstractVector{<:AbstractVector},symops::AbstractVector{<:AbstractSymop},maxiter = 100)
    positions = copy(sites)
    siteTypes = collect(eachindex(sites))
    i = 0
    for (x,site) in zip(siteTypes,positions)
        #new sites that are found are appended to positions and then used for generation iteratively. This loop should terminate eventually, since the number of sites in a unit cell is finite.
        #x is the type of the site, which does not change under symmetry transformations since, the site types are symmetry inequivalent by definition.
        for symop in symops
            r´ = symop(site) |> translate_back_to_UC
            if !approxin(r´,positions)
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

function getUnitCell(filename::AbstractString)
    sites = getproperty.(readVestaSites(filename),:position)
    symops = readVestaSymops(filename)
    getUnitCell(sites,symops)
end
"""given a filename return the Basis structure"""
function getBasis(filename::AbstractString)
    a1,a2,a3 = readLatticeVectors(filename)
    T = [a1 a2 a3]
    b,SiteType = getUnitCell(filename)
    
    b = [T*r for r in b]
    
    NUnique = maximum(SiteType)
    NNdist = min(getMinDistance(b),norm(a1),norm(a2),norm(a3))
    refSites_r = [T * r.position for r in readVestaSites(filename)]
    refSitePos = [findfirst(==(r), b) for r in refSites_r]
    refSites = Rvec.(0,0,0,refSitePos)
    return Basis_Struct_3D(;a1,a2,a3,b,NUnique,SiteType,NNdist,refSites)
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

"""Returns function acting on Rvec site by transforming and applying cartesian space function """
function gettransform(T::SiteTransformation,Basis::SpinFRGLattices.Basis_Struct)
    @inline transform(r::AbstractArray) = T(r)
    RV(x) = SpinFRGLattices.getRvec(x,Basis)
    @inline transform(R::Rvec) = getCartesian(R,Basis) |> T |> RV
    return transform
end

"""return the name of a file without the extension"""
function getName(filename)
    splitext(basename(filename))[1]
end

struct SiteTransformationRvec{T,B<:Basis_Struct} <: AbstractSymop
    T::SiteTransformation{T}
    Basis::B
end

using SpinFRGLattices:getRvec
(S::SiteTransformationRvec)(r::AbstractArray) = S.T(r)
(S::SiteTransformationRvec)(r::Rvec) = getRvec(S.T(getCartesian(r,S.Basis)),S.Basis)

function Base.show(io::IO, x::SiteTransformationRvec{T,B}) where {T,B}
    show(io,B)
    print(io," ")
    show(io,x.T)
end

Base.show(io::IO, ::MIME"text/plain", x::SiteTransformationRvec{T,B}) where {T,B} = show(io,x)

function splitSyms(syms::AbstractVector{T},refSites::AbstractArray) where T <: AbstractSymop
    
    nonRefSyms = T[]
    refSyms = [T[] for _ in eachindex(refSites)]

    for sym in syms
        isRefSym = false
        for (i,r) in enumerate(refSites)
            if sym(r) == r
                isRefSym = true
                push!(refSyms[i],sym)
            end
        end
        isRefSym || push!(nonRefSyms,sym)
    end
    return (;refSyms,nonRefSyms)
end

function getSymmetriesVesta(filename::AbstractString, Basis = getBasis(filename)::Basis_Struct)
    refSites = getCartesian.(Basis.refSites,Ref(Basis))
    
    syms = [SiteTransformationRvec(sym,Basis) for sym in readVestaSymops(filename)]
    return splitSyms(syms,refSites)
end


function reduceToInequivSites(NCell,nonRefSyms)
    InequivSites = Rvec_3D[]


    for b in 1:NCell
        R = Rvec(0,0,0,b)
        allRsymmetric = [translate_back_to_UC(sym(R)) for sym in nonRefSyms]

        if isempty(allRsymmetric ∩ InequivSites)
            push!(InequivSites,R)
        end
    end

    return InequivSites
end

import Base.*
*(S::SiteTransformation,S2::SiteTransformation) = SiteTransformation(S.WMatrix*S2.WMatrix)
*(S::SiteTransformationRvec,S2::SiteTransformationRvec) = SiteTransformationRvec(S.T*S2.T,S.Basis)

function readBondsVesta(filename::AbstractString)
    data = readFileInfo(filename,"SBOND",VestaFile())[begin:end-1]
    spl = getFields.(data)
    function parseline(s)
        _,site1,site2,minDistStr,maxDistStr,_ = s
        r,g,b = parse.(Int,s[end-2:end])
        colorRGB = (r,g,b)
        minDist = parse(Float64,minDistStr)
        maxDist = parse(Float64,maxDistStr)
        return (;site1,site2,minDist,maxDist,colorRGB)
    end
    bonds = parseline.(spl)
end

"""does not work for for lattices with more than one inequiv site yet, since the ref Symmetries need to be considered for each inequiv ref site individually"""
function generateSystem(NLen,filename;kwargs...)
    Name = getName(filename)
    Basis = getBasis(filename)
    (;refSyms,nonRefSyms) = getSymmetriesVesta(filename,Basis)
    System = getLatticeGeometry(NLen,Name,Basis,nonRefSyms,refSyms;kwargs...)
    return System
end

isInUnitCell(R::Rvec_3D) = R.n1 == R.n2 == R.n3 == 0

isidentity(s::SiteTransformation) = s.WMatrix == I
isidentity(s::SiteTransformationRvec) = s.T.WMatrix == I