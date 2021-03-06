# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    Rectangle{T} <: Surface{T}

Rectangular surface, not a valid CSG object.
The rotation of the rectangle around its normal is defined by `rotationvec`.
`rotationvecĂ—surfacenormal` is taken as the vector along the u axis.

**Can be used as a detector in [`AbstractOpticalSystem`](@ref)s.**

```julia
Rectangle(halfsizeu::T, halfsizev::T, [surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}]; rotationvec::SVector{3,T} = [0.0, 1.0, 0.0], interface::NullOrFresnel{T} = nullinterface(T))
```

The minimal case returns a rectangle centered at the origin with `surfacenormal = [0, 0, 1]`.
"""
struct Rectangle{T} <: PlanarShape{T}
    plane::Plane{T, 3}
    halfsizeu::T
    halfsizev::T
    uvec::SVector{3, T}
    vvec::SVector{3, T}

    function Rectangle(
        halfsizex::T,
        halfsizey::T;
        interface::NullOrFresnel{T} = NullInterface(T),
    ) where {T <: Real}
        @assert halfsizex > zero(T) && halfsizey > zero(T)
        return new{T}(
            Plane(
                zero(T),
                zero(T),
                one(T),
                zero(T),
                zero(T),
                zero(T),
                interface = interface,
            ),
            halfsizex,
            halfsizey,
            SVector{3, T}(1.0, 0.0, 0.0),
            SVector{3, T}(0.0, 1.0, 0.0),
        )
    end

    function Rectangle(
        halfsizeu::T,
        halfsizev::T,
        surfacenormal::AbstractArray{T, 1},
        centrepoint::AbstractArray{T, 1};
        interface::NullOrFresnel{T} = NullInterface(T),
        rotationvec::AbstractArray{T, 1} = SVector{3, T}(0.0, 1.0, 0.0),
    ) where {T <: Real}
        @assert length(surfacenormal) == 3 && length(centrepoint) == 3
        return Rectangle(
            halfsizeu,
            halfsizev,
            SVector{3, T}(surfacenormal),
            SVector{3, T}(centrepoint),
            interface = interface,
            rotationvec = SVector{3, T}(rotationvec),
        )
    end

    function Rectangle(
        halfsizeu::T,
        halfsizev::T,
        surfacenormal::SVector{3, T},
        centrepoint::SVector{3, T};
        interface::NullOrFresnel{T} = NullInterface(T),
        rotationvec::SVector{3, T} = SVector{3, T}(0.0, 1.0, 0.0),
    ) where {T <: Real}
        @assert halfsizeu > zero(T) && halfsizev > zero(T)
        nĚ‚ = normalize(surfacenormal)
        if abs(dot(rotationvec, nĚ‚)) == one(T)
            rotationvec = SVector{3, T}(1.0, 0.0, 0.0)
        end
        uvec = normalize(cross(normalize(rotationvec), nĚ‚))
        vvec = normalize(cross(nĚ‚, uvec))
        return new{T}(
            Plane(nĚ‚, centrepoint, interface = interface),
            halfsizeu,
            halfsizev,
            uvec,
            vvec,
        )
    end

    function Rectangle(
        plane::Plane{T, 3},
        halfsizeu::T,
        halfsizev::T,
        uvec::SVector{3, T},
        vvec::SVector{3, T},
    ) where {T <: Real}
        return new{T}(plane, halfsizeu, halfsizev, uvec, vvec)
    end
end

function Base.show(io::IO, a::Rectangle{T}) where {T <: Real}
    return print(
        io,
        "Rectangle{$T}($(centroid(a)), $(normal(a)), $(a.halfsizeu), $(a.halfsizev), $(interface(a)))",
    )
end

centroid(r::Rectangle{T}) where {T <: Real} = r.plane.pointonplane

function uvrange(::Type{Rectangle{T}}) where {T <: Real}
    return ((-one(T), one(T)), (-one(T), one(T)))
end

"""returns a 3D point. This takes into account the offset of centerpoint and the rotation vector used to construct the Rectangle. u and v are scaled by the size of the rectangle so that u=0,v=0 is one corner and u=v=1 is the diagonal corner. This function should go away once we have a sensible object transform hierarchy system."""
function point(r::Rectangle{T}, uvs::SMatrix{2, N, T}) where {N, T <: Real}
    result = MMatrix{3, N, T}(undef)
    for i in 1:N
        result[:, i] = point(r, uvs[1, i], uvs[2, i])
    end
    return SMatrix{3, N, T}(result)
end

"""returns a 3D point in the plane of the rectangle. This takes into account the offset of centerpoint and the rotation vector used to construct the Rectangle. u and v are scaled by the size of the rectangle so that u=0,v=0 is one corner and u=v=1 is the diagonal corner. This function should go away once we have a sensible object transform hierarchy system."""
function point(r::Rectangle{T}, u::T, v::T) where {T <: Real}
    return centroid(r) + (r.halfsizeu * u * r.uvec) + (r.halfsizev * v * r.vvec)
end
function partials(r::Rectangle{T}, ::T, ::T) where {T <: Real}
    return r.halfsizeu * r.uvec, r.halfsizev * r.vvec
end

function uv(r::Rectangle{T}, p::SVector{3, T}) where {T <: Real}
    return SVector{2, T}(
        dot(p - centroid(r), r.uvec) / r.halfsizeu,
        dot(p - centroid(r), r.vvec) / r.halfsizev,
    )
end

function onsurface(a::Rectangle{T}, point::SVector{3, T}) where {T <: Real}
    return onsurface(a.plane, point) &&
           abs(uv(a, point)[1]) <= one(T) &&
           abs(uv(a, point)[2]) <= one(T)
end

"""
    uvtopix(surf::Surface{T}, uv::SVector{2,T}, imsize::Tuple{Int,Int}) -> Tuple{Int,Int}

Converts a uvcoordinate on `surf` to an integer index to a pixel in an image of size `imsize`.
Not implemented on all `Surface` objects.
Used to determine where in the detector image a ray has hit when in intersects the detector surface of an [`AbstractOpticalSystem`](@ref).
"""
function uvtopix(
    ::Rectangle{T},
    uv::SVector{2, T},
    imsize::Tuple{Int, Int},
) where {T <: Real}
    u, v = uv
    h, w = imsize
    pixu = Int(floor((w - 1) * ((u + 1) / 2))) + 1
    pixv = h - Int(floor((h - 1) * ((v + 1) / 2)))
    return pixu, pixv
end

function surfaceintersection(
    rect::Rectangle{T},
    r::AbstractRay{T, 3},
) where {T <: Real}
    interval = surfaceintersection(rect.plane, r)
    if interval isa EmptyInterval{T} || isinfiniteinterval(interval)
        return EmptyInterval(T) # no ray plane intersection or inside plane but no hit
    else
        intersect = halfspaceintersection(interval)
        p = point(intersect)
        if abs(dot(p - centroid(rect), rect.uvec)) > rect.halfsizeu ||
           abs(dot(p - centroid(rect), rect.vvec)) > rect.halfsizev
            return EmptyInterval(T) # point outside rect
        else
            u, v = uv(rect, p)
            intuv = Intersection(
                Î±(intersect),
                p,
                normal(rect),
                u,
                v,
                interface(rect),
            )
            if dot(normal(rect), direction(r)) < zero(T)
                return positivehalfspace(intuv)
            else
                return rayorigininterval(intuv)
            end
        end
    end
end

"""returns the 2D vertices in the plane of the rectangle"""
function vertices(r::Rectangle{T}, ::Int = 0) where {T <: Real}
    return SMatrix{2, 4}(vertices3d(r)[1:2, :])
end

"""returns the vertices of the rectangle in 3D"""
function vertices3d(r::Rectangle{T}, ::Int = 0) where {T <: Real}
    pts = SVector{4, SVector{3, T}}(
        point(r, -one(T), -one(T)),
        point(r, -one(T), one(T)),
        point(r, one(T), one(T)),
        point(r, one(T), -one(T)),
    )
    temp = MMatrix{3, 4, T}(undef)

    for (j, pt) in pairs(pts)
        for i in 1:3
            temp[i, j] = pts[j][i]
        end
    end
    return SMatrix{3, 4, T}(temp)
end

function makemesh(r::Rectangle{T}, ::Int = 0) where {T <: Real}
    # p00,p01,p10,p11 = vertices(r)
    p00 = point(r, -one(T), -one(T))
    p01 = point(r, -one(T), one(T))
    p10 = point(r, one(T), -one(T))
    p11 = point(r, one(T), one(T))
    if validtri(p00, p11, p01) && validtri(p00, p10, p11)
        return TriangleMesh([Triangle(p00, p11, p01), Triangle(p00, p10, p11)])
    else
        return nothing
    end
end
