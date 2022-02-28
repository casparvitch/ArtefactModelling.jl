# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Surfaces

using FromFile: @from
using StaticArrays: SVector, SMatrix, MVector
using LinearAlgebra: cross, normalize, norm, dot

@from "Rays.jl" import Rays: AbstractRay, direction, origin, α
@from "Transforms.jl" import Transforms: Transform, identitytransform
@from "OpticalInterfaces.jl" import OpticalInterfaces:
    OpticalInterface,
    InterfaceMode,
    Reflect,
    Transmit,
    ReflectOrTransmit,
    NullInterface,
    FresnelInterface,
    opaqueinterface,
    NullOrFresnel,
    AllOpticalInterfaces

@from "BoundingBoxes.jl" import BoundingBoxes: BoundingBoxes, BoundingBox
@from "../utilities.jl" import samepoint,
    quadraticroots, NaNsafeatan, NaNsafeasin
@from "Intersections.jl" import Intersections:
    EmptyInterval,
    Intersection,
    rayorigininterval,
    Infinity,
    RayOrigin,
    isinfiniteinterval,
    halfspaceintersection,
    Interval

# export Primitive,
#     Surface,
#     surfaceintersection,
#     normal,
#     interface,
#     ParametricSurface,
#     point,
#     partials,
#     uvrange,
#     inside,
#     onsurface,
#     uv,
#     samplesurface,
#     triangulate,
#     makemesh

# NOTE All of the Primitives are include()'d at the bottom of this module

# from Plane.jl
# export Plane

# from Sphere.jl
# export Sphere

# from Cylinder.jl
# export Cylinder

# from SphericalCap.jl
# export SphericalCap

# from NonCSGs
# export Triangle,
#     TriangleMesh,
#     PlanarShape,
#     Rectangle,
#     vertices3d,
#     Hexagon,
#     ConvexPolygon,
#     localframe,
#     Ellipse,
#     Circle,
#     InfiniteStop,
#     FiniteStop,
#     RectangularAperture,
#     CircularAperture,
#     Annulus,
#     InfiniteStopConvexPoly,

# from Curves
# export BezierCurve,
#     BezierSurface,
#     BSplineCurve,
#     BSplineSurface,
#     KnotVector,
#     PowerBasisCurve,
#     CurveType,
#     Spline,
#     SplineSurface

# from AsphericSurface.jl
# export AsphericSurface,
#     EvenAsphericSurface,
#     OddAsphericSurface,
#     OddEvenAsphericSurface,
#     asphericType,
#     EVEN,
#     CONIC,
#     ODD,
#     ODDEVEN

# from AccelSurface.jl
# export AcceleratedParametricSurface

# from Primitives/{Zernike,Qtype,Chebyshev}.jl
# export ZernikeSurface,
#     ZernikeIndexType,
#     ZernikeIndexingOSA,
#     ZernikeIndexingNoll,
#     QTypeSurface,
#     ChebyshevSurface,

"""
    Primitive{T<:Real}

`T` is the number type used to represent the primitive,  e.g., `Float64`.
Primitives are the basic elements which can be stored in bounding volume hierarchies and include surfaces and CSG objects

**Must** implement the following:
```julia
boundingbox(a::Primitive{T})::BoundingBox{T}
centroid(a::Primitive{T})::SVector{3,T}
```
"""
abstract type Primitive{T <: Real} end

"""
    Surface{T<:Real}

`T` is the number type used to represent the surface, e.g., `Float64`.
Basic `Surface`s are _not_ valid CSG objects, they function only in a stand-alone capacity.

**Must** implement the following:
```julia
surfaceintersection(surface::Surface{T}, ray::AbstractRay{T,3}) -> Union{EmptyInterval{T},Interval{T}}
normal(surface::Surface{T}) -> SVector{3,T}
interface(surface::Surface{T}) -> OpticalInterface{T}
makemesh(surface::Surface{T}) -> TriangleMesh{T}
```

In a conventional ray tracer the surface intersection function would only return the first surface the ray intersects. Because our ray tracer does CSG operations the surface intersection function intersects the ray with all leaf surfaces which are part of the CSG tree. 

Each leaf surface returns one or more 1D intervals along the ray. These intervals contain the part of the ray which is inside the surface. The intervals computed at the leaves are propagated upward through the CSG tree and the CSG operations of union, intersection, and difference are applied to generate new intervals which are themselves propagated upward.

The result is a union of 1D intervals, which may be disjoint, a single interval, or empty. The union of intervals represents the parts of the ray which are inside the CSG object.

Inside is well defined for halfspaces such as cylinders and spheres which divide space into two parts, but not for Bezier or NURBS patches which generally do not enclose a volume.  For surfaces which are not halfspaces the notion of inside is defined locally by computing the angle between the incoming ray and the normal of the surface at the point of intersection. All surfaces must be defined so that the normal points to the outside of the surface. 

A negative dot product between the incoming ray and the normal indicates the ray is coming from the outside of the surface and heading toward the inside. A positive dot product indicates the ray is coming from the inside of the surface and heading toward the outside.

Intervals are defined along the ray which is being intersected with the surface, so they are one dimensional. For example, assume we have a ray with origin o on the outside of a plane and an intersection with the plane at point int = o + td where t is a scalar and d is the unit direction of the ray. The inside interval will be (Intersection(t),Infinity). This interval begins at the intersection point on the plane and continues to positive infinity. The Intersection struct stores both the parametric value t and the 3D point of intersection to make various operations more efficient. But the interval operations only depend on the parametric value t.

If the origin o is on the inside of the plane then the inside interval will be (RayOrigin,Intersection(t)). Only the part of the ray from the ray origin to the intersection point is inside the plane. 

It is the programmer's responsibility to return Interval results from surfaceintersection that maintain these properties.

The following must be impemented only if the surface is being used as a detector
```julia
uv(surface::Surface{T}, p::SVector{3,T}) -> SVector{2,T}
uvtopix(surface::Surface{T}, uv::SVector{2,T}, imsize::Tuple{Int,Int}) -> Tuple{Int,Int}
onsurface(surface::Surface{T}, p::SVector{3,T}) -> Bool
```
"""
abstract type Surface{T <: Real} <: Primitive{T} end

"""
    ParametricSurface{T,N} <: Surface{T}

`T` is the number type used to represent the surface, e.g., `Float64`.
`N` is the dimension of the space the surface is embedded in.
`ParametricSurface`s are valid CSG objects, in some cases (where analytic intersection isn't possible) they must be wrapped in an [`AcceleratedParametricSurface`](@ref) for use.

**Must** implement the following:
```julia
uv(surface::ParametricSurface{T,N}, p::SVector{N,T}) -> SVector{2,T}
uvrange(surface::ParametricSurface{T,N}) -> Tuple{Tuple{T,T},Tuple{T,T}}
point(surface::ParametricSurface{T,N}, u::T, v::T) -> SVector{N,T}
partials(surface::ParametricSurface{T,N}, u::T, v::T) -> Tuple{SVector{N,T}, SVector{N,T}}
normal(surface::ParametricSurface{T,N}, u::T, v::T) -> SVector{N,T}
inside(surface::ParametricSurface{T,N}, p: :SVector{N,T}) -> Bool
onsurface(surface::ParametricSurface{T,N}, p::SVector{N,T}) -> Bool
surfaceintersection(surface::ParametricSurface{T,N}, AbstractRay::Ray{T,N}) -> Union{EmptyInterval{T},Interval{T},DisjointUnion{T}}
interface(surface::ParametricSurface{T,N}) -> OpticalInterface{T}
```
"""
abstract type ParametricSurface{S <: Real, N} <: Surface{S} end

# chuck this in here to maintain a nice heirarchy
function BoundingBoxes.BoundingBox(
    s::ParametricSurface{T, 3},
    transform::Transform{T} = identitytransform(T),
) where {T <: Real}
    # get the bounding box of a transformed bounding box
    bbox = BoundingBox(s)
    if transform == identitytransform(T)
        return bbox
    else
        p1 = transform * SVector(bbox.xmin, bbox.ymin, bbox.zmin)
        p2 = transform * SVector(bbox.xmin, bbox.ymax, bbox.zmin)
        p3 = transform * SVector(bbox.xmin, bbox.ymax, bbox.zmax)
        p4 = transform * SVector(bbox.xmin, bbox.ymin, bbox.zmax)
        p5 = transform * SVector(bbox.xmax, bbox.ymin, bbox.zmin)
        p6 = transform * SVector(bbox.xmax, bbox.ymax, bbox.zmin)
        p7 = transform * SVector(bbox.xmax, bbox.ymax, bbox.zmax)
        p8 = transform * SVector(bbox.xmax, bbox.ymin, bbox.zmax)
        return BoundingBox(SVector(p1, p2, p3, p4, p5, p6, p7, p8))
    end
end

# all subclasses must implement one of these at least...
"""
    point(surf::ParametricSurface{T}, u::T, v::T) -> SVector{3,T}
    point(surf::ParametricSurface{T}, uv::SVector{2,T}) -> SVector{3,T}

Returns the 3D point on `surf` at the given uv coordinate.
"""
function point(s::ParametricSurface{T}, u::T, v::T) where {T <: Real}
    return point(s, SVector{2, T}(u, v))
end
function point(s::ParametricSurface{T}, uv::SVector{2, T}) where {T <: Real}
    return point(s, uv[1], uv[2])
end

"""
    point(ray::AbstractRay{T,N}, alpha::T) -> SVector{T, N}

Returns a point on the ray at origin + alpha * direction. Alpha must be >= 0.
"""
function point(ray::AbstractRay{T, N}, alpha::T) where {N, T <: Real}
    @assert alpha >= zero(T) "Alpha must be nonnegative. alpha value: $alpha"
    return origin(ray) + alpha * direction(ray)
end
point(a::Intersection{T, N}) where {T <: Real, N} = a.point

"""
    normal(surf::ParametricSurface{T}, u::T, v::T) -> SVector{3,T}
    normal(surf::ParametricSurface{T}, uv::SVector{2,T}) -> SVector{3,T}

Returns the normal to `surf` at the given uv coordinate.
"""
function normal(s::ParametricSurface{T}, u::T, v::T) where {T <: Real}
    return normal(s, SVector{2, T}(u, v))
end
function normal(s::ParametricSurface{T}, uv::SVector{2, T}) where {T <: Real}
    return normal(s, uv[1], uv[2])
end
"""
    partials(surf::ParametricSurface{T}, u::T, v::T) -> (SVector{3,T}, SVector{3,T})
    partials(surf::ParametricSurface{T}, uv::SVector{2,T}) -> (SVector{3,T}, SVector{3,T})

Returns a tuple of the 3D partial derivatives of `surf` with respect to u and v at the given uv coordinate.
"""
function partials(s::ParametricSurface{T}, u::T, v::T) where {T <: Real}
    return partials(s, SVector{2, T}(u, v))
end
function partials(s::ParametricSurface{T}, uv::SVector{2, T}) where {T <: Real}
    return partials(s, uv[1], uv[2])
end
"""
    uv(surf::ParametricSurface{T}, p::SVector{3,T}) -> SVector{2,T}
    uv(surf::ParametricSurface{T}, x::T, y::T, z::T) -> SVector{2,T}

Returns the uv coordinate on `surf` of a point, `p`, in 3D space.
If `onsurface(surf, p)` is false then the behavior is undefined, it may return an inorrect uv, an invalid uv, NaN or crash.
"""
function uv(s::ParametricSurface{T, 3}, x::T, y::T, z::T) where {T <: Real}
    return uv(s, SVector{3, T}(x, y, z))
end
function uv(s::ParametricSurface{T, 3}, p::SVector{3, T}) where {T <: Real}
    return uv(s, p[1], p[2], p[3])
end
"""
    inside(surf::ParametricSurface{T}, p::SVector{3,T}) -> Bool
    inside(surf::ParametricSurface{T}, x::T, y::T, z::T) -> Bool

Tests whether a 3D point in world space is _inside_ `surf`.
"""
function inside(s::ParametricSurface{T, 3}, x::T, y::T, z::T) where {T <: Real}
    return inside(s, SVector{3, T}(x, y, z))
end
function inside(s::ParametricSurface{T, 3}, p::SVector{3, T}) where {T <: Real}
    return inside(s, p[1], p[2], p[3])
end
"""
    onsurface(surf::ParametricSurface{T}, p::SVector{3,T}) -> Bool
    onsurface(surf::ParametricSurface{T}, x::T, y::T, z::T) -> Bool

Tests whether a 3D point in world space is _on_ `surf`.
"""
function onsurface(
    s::ParametricSurface{T, 3},
    x::T,
    y::T,
    z::T,
) where {T <: Real}
    return onsurface(s, SVector{3, T}(x, y, z))
end
function onsurface(
    s::ParametricSurface{T, 3},
    p::SVector{3, T},
) where {T <: Real}
    return onsurface(s, p[1], p[2], p[3])
end
"""
    uvrange(s::ParametricSurface)
    uvrange(::Type{S}) where {S<:ParametricSurface}

Returns a tuple of the form: `((umin, umax), (vmin, vmax))` specifying the limits of the parameterisation for this surface type.
Also implemented for some `Surface`s which are not `ParametricSurface`s (e.g. `Rectangle`).
"""
uvrange(::S) where {S <: ParametricSurface} = uvrange(S)

"""
    samplesurface(surf::ParametricSurface{T,N}, samplefunction::Function, numsamples::Int = 30)

Sample a parametric surface on an even `numsamples`×`numsamples` grid in UV space with provided function
"""
function samplesurface(
    surf::ParametricSurface{T, N},
    samplefunction::Function,
    numsamples::Int = 30,
) where {T <: Real, N}
    urange, vrange = uvrange(surf)
    rt = typeof(samplefunction(surf, urange[1], vrange[1]))
    samples = Vector{rt}(undef, (numsamples + 1)^2)
    ustep = (urange[2] - urange[1]) / numsamples
    vstep = (vrange[2] - vrange[1]) / numsamples
    for ui in 0:numsamples
        for vi in 0:numsamples
            u = urange[1] + ui * ustep
            v = vrange[1] + vi * vstep
            samples[ui * (numsamples + 1) + vi + 1] = samplefunction(surf, u, v)
        end
    end
    return samples
end

# order is important here (for some)
include("Primitives/Plane.jl")
include("Primitives/Sphere.jl")
include("Primitives/Cylinder.jl")
include("Primitives/SphericalCap.jl")

include("Primitives/NonCSG/Triangle.jl")
include("Primitives/NonCSG/PlanarShape.jl")
include("Primitives/NonCSG/Rectangle.jl")
include("Primitives/NonCSG/Hexagon.jl")
include("Primitives/NonCSG/ConvexPolygon.jl")
include("Primitives/NonCSG/Ellipse.jl")
include("Primitives/NonCSG/Stop.jl")

include("AccelSurface.jl")
include("Primitives/AsphericSurface.jl")
include("Primitives/Zernike.jl")
include("Primitives/Qtype.jl")
include("Primitives/Chebyshev.jl")

include("Primitives/Curves/Knots.jl")
include("Primitives/Curves/Spline.jl")
include("Primitives/Curves/BSpline.jl")
include("Primitives/Curves/Bezier.jl")
include("Primitives/Curves/PowerBasis.jl")

end # module
