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
function readUniqueSites(filename)
    Info = readFileInfo(filename,"STRUC",VestaFile())[1:2:end-1] # skip standard uncertainty lines and line of trailing zeros
    function parseSites(s)
        number,atom,label,occ,x,y,z,_ = getFields(s)
        position = parse.(Float64,SA[x,y,z])
        return (;atom,label,position)
    end
    return parseSites.(Info)
end


abstract type AbstractSymop <: Function end

struct FractionalTransformation{T<:Real} <: AbstractSymop
    WMatrix::SMatrix{4,4,T,16}
end

getTranslation(S::FractionalTransformation) = S.WMatrix[1:3,4]
getRotation(S::FractionalTransformation) = SMatrix{3,3}(S.WMatrix[1:3,1:3])


FractionalTransformation(t::SVector{3,T},W::SMatrix{3,3,T,9}) where T <: Real = FractionalTransformation(SA[
    W[1,1] W[1,2] W[1,3] t[1];
    W[2,1] W[2,2] W[2,3] t[2];
    W[3,1] W[3,2] W[3,3] t[3];
    0 0 0 1
    ]
) 

import Base.*
*(S::FractionalTransformation,S2::FractionalTransformation) = FractionalTransformation(S.WMatrix*S2.WMatrix)

import Base.^
^(S::FractionalTransformation,i::Number) = FractionalTransformation(S.WMatrix^i)

"""given an origin and a transformation matrix, returns a FractionalTransformation"""
function fractionaltransformation_origin(origin::AbstractVector{T},matrix::AbstractMatrix{T}) where T <: Real
    return FractionalTransformation(-matrix*origin +origin,matrix)
end

# getOrigin(S::FractionalTransformation) = inv(getRotation(S)+I)*getTranslation(S)

function (S::FractionalTransformation)(vec::AbstractVector{T}) where {T<:Real}
    x,y,z = vec
    x,y,z,_ = S.WMatrix * SA[x,y,z,one(x)]
    return SA[x,y,z]
end

parseAsSMatrix(v::AbstractVector{<:AbstractString}) = SMatrix{3,3,Float64,9}([parse(Float64,s) for s in v])'

function Base.show(io::IO, x::FractionalTransformation{T}) where T
    println(io,"3 dim SiteTransformation{",T,"}")
    Base.print_matrix(io,x.WMatrix)
end
Base.show(io::IO, ::MIME"text/plain", x::FractionalTransformation) = show(io,x)

ispuretranslation(S::FractionalTransformation) = abs(det(getRotation(S))) == 0

"""read the symmetry operations from a vesta file"""
function readSymops(filename)
    Info = readFileInfo(filename,"SYMOP",VestaFile())
    Info = getFields.(Info)
    origins = [parseAsSVector(data[1:3]) for data in Info]
    matrices = [parseAsSMatrix(data[4:12]) for data in Info]
    Syms = FractionalTransformation.(origins,matrices)
    filter!(!ispuretranslation,Syms)
    return Syms
end

"""for a position in fractional coordinates return true if it is in the unit cell"""
function isInUnitCell(pos_fractional::AbstractVector,lower=0,upper=1)
    all(pos_fractional .>= lower) && all(pos_fractional .< upper)
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
    sites = getproperty.(readUniqueSites(filename),:position)
    symops = readSymops(filename)
    getUnitCell(sites,symops)
end
"""given a filename return the Basis structure"""
function getBasis(filename::AbstractString)
    a1,a2,a3 = readLatticeVectors(filename)
    T = [a1 a2 a3]
    b,SiteType = getUnitCell(filename)
    
    b = [T*r for r in b]
    
    NUnique = maximum(SiteType)
    NNdist = min(getMinDistance(b),norm(a1),norm(a2),norm(a3))# not entirely correct, need to generate more than one unit cell
    refSites_r = [T * r.position for r in readUniqueSites(filename)]
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

"""return the name of a file without the extension"""
function getName(filename)
    splitext(basename(filename))[1]
end

struct SiteTransformation{Ty,B<:Basis_Struct} <: AbstractSymop
    T::FractionalTransformation{Ty}
    Basis::B
end

import Base.*
*(S::SiteTransformation,S2::SiteTransformation) = SiteTransformation(S.T*S2.T,S.Basis)

import Base.^
^(S::SiteTransformation,i::Number) = SiteTransformation(S.T^i,S.Basis)

"""returns a SiteTransformation accounting for the correct coordinates of the basis"""
function siteTransformation(S::FractionalTransformation{Ty},Basis::B) where {Ty,B<:Basis_Struct}
    T = Basis.T
    Tinv = inv(T)
    transl = getTranslation(S)
    mat = getRotation(S)
    s = FractionalTransformation(Tinv*transl,Tinv*mat*T)
    return SiteTransformation(s,Basis)
end
using SpinFRGLattices:getRvec
(S::SiteTransformation)(r::AbstractArray) = S.T(r)
(S::SiteTransformation)(r::Rvec) = getRvec(S.T(getCartesian(r,S.Basis)),S.Basis)

function Base.show(io::IO, x::SiteTransformation{T,B}) where {T,B}
    show(io,B)
    print(io," ")
    show(io,x.T)
end

Base.show(io::IO, ::MIME"text/plain", x::SiteTransformation{T,B}) where {T,B} = show(io,x)

function splitSyms(syms::AbstractVector{T},refSites::AbstractArray) where T <: AbstractSymop
    
    nonRefSyms = T[]
    refSyms = [T[] for _ in eachindex(refSites)]

    for sym in syms
        isRefSym = false
        for (i,r) in enumerate(refSites)
            if sym(r) ≈ r
                isRefSym = true
                push!(refSyms[i],sym)
            end
        end
        isRefSym || push!(nonRefSyms,sym)
    end
    return (;refSyms,nonRefSyms)
end

getTranslation(S::SiteTransformation{Float64, Basis_Struct_3D}) = getTranslation(S.T)

function getSymmetries(filename::AbstractString, Basis = getBasis(filename)::Basis_Struct)
    refSites = [getCartesian(r,Basis) for r in Basis.refSites]
    syms = [siteTransformation(sym,Basis) for sym in readSymops(filename)]
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

function readBonds(filename::AbstractString)
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

function generateSiteSyms(syms::AbstractVector{<:FractionalTransformation},Basis)
    Generatedsyms = generateSymms(syms)
    # Generatedsyms = syms
    Rsyms = [siteTransformation(sym,Basis) for sym in Generatedsyms]
    refSites = [getCartesian(r,Basis) for r in Basis.refSites]
    return splitSyms(Rsyms,refSites)
end

"""does not work for for lattices with more than one inequiv site yet, since the ref Symmetries need to be considered for each inequiv ref site individually"""
function generateSystem(NLen,filename;addSyms = nothing,kwargs...)
    Name = getName(filename)*"_NLen=$(NLen)"
    Basis = getBasis(filename)

    (;refSyms,nonRefSyms) = generateSymmetries(filename,addSyms;Basis)
    System = getLatticeGeometry(NLen,Name,Basis,nonRefSyms,refSyms;kwargs...)
    return System
end

function generateSymmetries(filename::AbstractString,addSyms=nothing; Basis = getBasis(filename))
    syms = readSymops(filename)

    if addSyms !== nothing
        append!(syms,addSyms)
    end
    return generateSiteSyms(syms,Basis)
        
end

isInUnitCell(R::Rvec_3D) = R.n1 == R.n2 == R.n3 == 0

isidentity(s::FractionalTransformation) = s.WMatrix == I
isidentity(s::SiteTransformation) = s.T.WMatrix == I

import SpinFRGLattices.generateReducedLattice

function generateReducedLattice(NLen,filename,addSyms=nothing,kwargs...)
    Basis = getBasis(filename)
    (;refSyms,nonRefSyms) = generateSymmetries(filename,addSyms;Basis)

    return generateReducedLattice(NLen,Basis,nonRefSyms,refSyms;kwargs...)

end
