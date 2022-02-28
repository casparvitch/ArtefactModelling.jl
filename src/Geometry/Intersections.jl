
module Intersections

using FromFile: @from
using StaticArrays: SVector, SMatrix
using LinearAlgebra: cross, normalize, norm, dot
@from "OpticalInterfaces.jl" import OpticalInterfaces:
    #     OpticalInterface,
    #     InterfaceMode,
    #     Reflect,
    #     Transmit,
    #     ReflectOrTransmit,
    #     NullInterface,
    #     FresnelInterface,
    #     opaqueinterface,
    #     NullOrFresnel,
    AllOpticalInterfaces
@from "Transforms.jl" import Transforms: Transform, identitytransform, rotate
@from "Rays.jl" import Rays

# ============================================================================ #
# 88 88b 88 888888 888888 88""Yb .dP"Y8 888888  dP""b8 888888 88  dP"Yb  88b 88 
# 88 88Yb88   88   88__   88__dP `Ybo." 88__   dP   `"   88   88 dP   Yb 88Yb88 
# 88 88 Y88   88   88""   88"Yb  o.`Y8b 88""   Yb        88   88 Yb   dP 88 Y88 
# 88 88  Y8   88   888888 88  Yb 8bodP' 888888  YboodP   88   88  YbodP  88  Y8 
# ============================================================================ #

"""
Each [`Interval`](@ref) consists of two `IntervalPoint`s, one of 
[`RayOrigin`](@ref), [`Intersection`](@ref) or [`Infinity`](@ref).
"""
abstract type IntervalPoint{T <: Real} end
Base.eltype(::IntervalPoint{T}) where {T <: Real} = T

abstract type FinitePoint{T} <: IntervalPoint{T} end
Base.eltype(::FinitePoint{T}) where {T <: Real} = T

"""
    Intersection{T,N} <: IntervalPoint{T}

Represents the point at which an [`Ray`](@ref) hits a [`Surface`](@ref).
This consists of the distance along the ray, the intersection point in world
space, the normal in world space, the UV on the surface and the [`OpticalInterface`](@ref) hit.

Has the following accessor methods:
```julia
point(a::Intersection{T,N}) -> SVector{N,T} # NOTE: this is in Surfaces.jl now
normal(a::Intersection{T,N}) -> SVector{N,T}
uv(a::Intersection{T,N}) -> SVector{2,T}
u(a::Intersection{T,N}) -> T
v(a::Intersection{T,N}) -> T
α(a::Intersection{T,N}) -> T
interface(a::Intersection{T,N}) -> OpticalInterface{T}
flippednormal(a::Intersection{T,N}) -> Bool
```
"""
struct Intersection{T, N} <: FinitePoint{T}
    α::T # value of ray parameter at the point of intersection
    point::SVector{N, T}
    normal::SVector{N, T}
    u::T
    v::T
    interface::AllOpticalInterfaces{T} # returns a union of all OpticalInterface
    # subtypes - can't have an abstract type here as it results in allocations
    flippednormal::Bool

    function Intersection(
        α::T,
        point::SVector{N, T},
        normal::SVector{N, T},
        u::T,
        v::T,
        interface::AllOpticalInterfaces{T};
        flippednormal = false,
    ) where {T <: Real, N}
        return new{T, N}(
            α,
            point,
            normalize(normal),
            u,
            v,
            interface,
            flippednormal,
        )
    end

    function Intersection(
        α::T,
        point::AbstractVector{T},
        normal::AbstractVector{T},
        u::T,
        v::T,
        interface::AllOpticalInterfaces{T};
        flippednormal = false,
    ) where {T <: Real}
        @assert length(point) == length(normal)
        N = length(point)
        return new{T, N}(
            α,
            SVector{N, T}(point),
            SVector{N, T}(normal),
            u,
            v,
            interface,
            flippednormal,
        )
    end
end

normal(a::Intersection{T, N}) where {T <: Real, N} = a.normal
uv(a::Intersection{T, N}) where {T <: Real, N} = SVector{2, T}(a.u, a.v)
u(a::Intersection{T, N}) where {T <: Real, N} = a.u
v(a::Intersection{T, N}) where {T <: Real, N} = a.v
Rays.α(a::Intersection{T, N}) where {T <: Real, N} = a.α
interface(a::Intersection{T, N}) where {T <: Real, N} = a.interface
flippednormal(a::Intersection{T, N}) where {T <: Real, N} = a.flippednormal

function Base.print(io::IO, a::Intersection{T, N}) where {T <: Real, N}
    println(io, Intersection{T, N})
    println(io, "α \t$(Rays.α(a))")
    println(io, "Point \t$(a.point)")
    println(io, "Normal \t$(normal(a))")
    println(io, "normal is flipped? \t$(flippednormal(a))")
    println(io, "u \t$(u(a))")
    println(io, "v \t$(v(a))")
    return println(io, "interface \t$(interface(a))")
end

"""
    reversenormal(a::Intersection{T,N})

Used by the CSG complement operator (i.e. [`-`](@ref)) to reverse the inside 
outside sense of the object.
"""
function reversenormal(a::Intersection{T, N}) where {T <: Real, N}
    return Intersection(
        Rays.α(a),
        a.point,
        -normal(a),
        u(a),
        v(a),
        interface(a),
        flippednormal = !flippednormal(a),
    )
end

"""
    Infinity{T} <: IntervalPoint{T}

Point representing ∞ within an [`Interval`](@ref).

```julia
Infinity(T = Float64)
Infinity{T}()
```
"""
struct Infinity{T} <: IntervalPoint{T}
    Infinity(::Type{T} = Float64) where {T <: Real} = new{T}()
    Infinity{T}() where {T <: Real} = new{T}()
end
"""
    RayOrigin{T} <: IntervalPoint{T}

Point representing 0 within an [`Interval`](@ref), i.e. the start of the ray.

```julia
RayOrigin(T = Float64)
RayOrigin{T}()
```
"""
struct RayOrigin{T} <: FinitePoint{T}
    RayOrigin(::Type{T} = Float64) where {T <: Real} = new{T}()
    RayOrigin{T}() where {T <: Real} = new{T}()
end

Base.:(==)(::RayOrigin{P}, ::RayOrigin{P}) where {P <: Real} = true
Base.:(==)(::Infinity{P}, ::Infinity{P}) where {P <: Real} = true
function Base.:(==)(
    a::Intersection{P, N},
    b::Intersection{P, N},
) where {P <: Real, N}
    return Rays.α(a) == Rays.α(b)
end
Base.:(==)(::IntervalPoint{P}, ::IntervalPoint{P}) where {P <: Real} = false

Base.:(<=)(::RayOrigin{P}, ::RayOrigin{P}) where {P <: Real} = true
Base.:(<=)(::Infinity{P}, ::Infinity{P}) where {P <: Real} = true
Base.:(<=)(a::IntervalPoint{P}, b::Infinity{P}) where {P <: Real} = true
function Base.:(<=)(
    a::Intersection{P, N},
    b::Intersection{P, N},
) where {P <: Real, N}
    return Rays.α(a) <= Rays.α(b)
end
Base.:(<=)(a::RayOrigin{P}, b::Infinity{P}) where {P <: Real} = true
Base.:(<=)(a::RayOrigin{P}, b::IntervalPoint{P}) where {P <: Real} = true
Base.:(<=)(::IntervalPoint{P}, ::IntervalPoint{P}) where {P <: Real} = false

Base.:(>=)(::RayOrigin{P}, ::RayOrigin{P}) where {P <: Real} = true
Base.:(>=)(::Infinity{P}, ::Infinity{P}) where {P <: Real} = true
Base.:(>=)(a::IntervalPoint{P}, b::Infinity{P}) where {P <: Real} = false
function Base.:(>=)(
    a::Intersection{P, N},
    b::Intersection{P, N},
) where {P <: Real, N}
    return Rays.α(a) >= Rays.α(b)
end
Base.:(>=)(a::RayOrigin{P}, b::Infinity{P}) where {P <: Real} = false
Base.:(>=)(a::RayOrigin{P}, b::IntervalPoint{P}) where {P <: Real} = false
Base.:(>=)(::IntervalPoint{P}, ::IntervalPoint{P}) where {P <: Real} = false

Base.:(<)(::Infinity{P}, ::Infinity{P}) where {P <: Real} = false
Base.:(<)(::IntervalPoint{P}, ::Infinity{P}) where {P <: Real} = true
function Base.:(<)(
    a::Intersection{P, N},
    b::Intersection{P, N},
) where {P <: Real, N}
    return Rays.α(a) < Rays.α(b)
end

Base.:(>)(::Infinity{P}, ::Infinity{P}) where {P <: Real} = false
Base.:(>)(::IntervalPoint{P}, ::Infinity{P}) where {P <: Real} = false
function Base.:(>)(
    a::Intersection{P, N},
    b::Intersection{P, N},
) where {P <: Real, N}
    return Rays.α(a) > Rays.α(b)
end

function Base.isless(a::IntervalPoint{P}, b::IntervalPoint{P}) where {P <: Real}
    return Rays.α(a) < Rays.α(b)
end

Base.eltype(::Intersection{T, N}) where {T <: Real, N} = T

"""
    isinfinity(a) -> Bool

Returns true if `a` is [`Infinity`](@ref). In performance critical contexts use 
`a isa Infinity{T}`.

"""
isinfinity(::Infinity) = true
isinfinity(::Any) = false

Rays.α(::Infinity{T}) where {T <: AbstractFloat} = typemax(T)
Base.eltype(::Infinity{T}) where {T <: Real} = T

"""
    israyorigin(a) -> Bool

Returns true if `a` is [`RayOrigin`](@ref). In performance critical contexts use 
`a isa RayOrigin{T}`.

"""
israyorigin(::Any) = false
israyorigin(::RayOrigin) = true

Rays.α(::RayOrigin{T}) where {T <: Real} = zero(T)
Base.eltype(::RayOrigin{T}) where {T <: Real} = T

"""
Apply a Transform to an Intersection object
"""
function Base.:*(
    a::Transform{T},
    intsct::Intersection{T, 3},
)::Intersection{T, 3} where {T <: Real}
    u, v = uv(intsct)
    i = interface(intsct)
    return Intersection(
        Rays.α(intsct),
        a * intsct.point,
        rotate(a, normal(intsct)),
        u,
        v,
        interface(intsct),
        flippednormal = flippednormal(intsct),
    )
end

# ============================================================================ #
#        88 88b 88 888888 888888 88""Yb Yb    dP    db    88     
#        88 88Yb88   88   88__   88__dP  Yb  dP    dPYb   88     
#        88 88 Y88   88   88""   88"Yb    YbdP    dP__Yb  88  .o 
#        88 88  Y8   88   888888 88  Yb    YP    dP""""Yb 88ood8 
# ============================================================================ #

# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

abstract type AbstractRayInterval{T <: Real} end

"""
    EmptyInterval{T} <: AbstractRayInterval{T}

An interval with no [`Intersection`](@ref)s which is also not infinite.

```
EmptyInterval(T = Float64)
EmptyInterval{T}()
```
"""
struct EmptyInterval{T} <: AbstractRayInterval{T}
    EmptyInterval(::Type{T} = Float64) where {T <: Real} = new{T}()
    EmptyInterval{T}() where {T <: Real} = new{T}()
end

"""
    Interval{T} <: AbstractRayInterval{T}

Datatype representing an interval between two [`IntervalPoint`](@ref)s on a ray.

The lower element can either be [`RayOrigin`](@ref) or an [`Intersection`](@ref).
The upper element can either be an [`Intersection`](@ref) or [`Infinity`](@ref).

```julia
positivehalfspace(int::Intersection) -> 
    Interval with lower = int, upper = Infinity
rayorigininterval(int::Intersection) -> 
    Interval with lower = RayOrigin, upper = int
Interval(low, high)
```

Has the following accessor methods:
```julia
lower(a::Interval{T}) -> Union{RayOrigin{T},Intersection{T,3}}
upper(a::Interval{T}) -> Union{Intersection{T,3},Infinity{T}}
```
"""
struct Interval{R} <: AbstractRayInterval{R}
    lower::Union{RayOrigin{R}, Intersection{R, 3}} #making this a Union of two 
    # concrete types allows the compiler to statically enough space to hold the 
    # larger of the two. The entire Interval struct can then be stack allocated.
    upper::Union{Intersection{R, 3}, Infinity{R}} #same as above

    function Interval(
        low::S,
        high::T,
    ) where {
        R <: Real,
        S <: Union{RayOrigin{R}, Intersection{R, 3}},
        T <: Union{Intersection{R, 3}, Infinity{R}},
    }
        @assert low <= high
        return new{R}(low, high)
    end
end

function Base.show(io::IO, a::Interval{R}) where {R <: Real}
    println(io, Interval{R})
    println(io, "\tLower \n\t\t$(lower(a))")
    return println(io, "\tUpper \n\t\t$(upper(a))")
end

Base.eltype(::Interval{R}) where {R} = R

lower(a::Interval) = a.lower
upper(a::Interval) = a.upper

# ============================================================================ #

"""
To prevent allocations we have a manually managed pool of arrays of 
[`Interval`](@ref)s which are used to store values during execution.
The memory is kept allocated and reused across runs of functions like 
[`trace`](@ref).

`threadedintervalpool` is a global threadsafe pool which is accessed through 
the functions:
```julia
newinintervalpool!(::Type{T} = Float64, tid::Int = Threads.threadid(
    )) -> Vector{Interval{T}}
indexednewinintervalpool!(::Type{T} = Float64, tid::Int = Threads.threadid(
    )) -> Tuple{Int,Vector{Interval{T}}}
emptyintervalpool!(::Type{T} = Float64, tid::Int = Threads.threadid())
getfromintervalpool([::Type{T} = Float64], id::Int, tid::Int = Threads.threadid(
    )) -> Vector{Interval{T}}
```
"""
struct IntervalPool{T <: Real}
    allocated::Vector{Vector{Interval{T}}}
    unallocated::Vector{Vector{Interval{T}}}

    function IntervalPool{T}() where {T <: Real}
        return new{T}(
            Vector{Vector{Interval{T}}}(),
            Vector{Vector{Interval{T}}}(),
        )
    end
end

function allocate!(a::IntervalPool{T}) where {T <: Real}
    if length(a.unallocated) === 0
        #this will extend the array with a new Vector of intervals.
        push!(a.unallocated, Vector{Interval{T}}())
    end

    temp = pop!(a.unallocated)
    Base.empty!(temp) #resets count in array but doesn't reclaim space
    push!(a.allocated, temp)
    return temp
end

function empty!(a::IntervalPool{T}) where {T <: Real}
    while length(a.allocated) > 0
        temp = pop!(a.allocated)
        Base.empty!(temp) #not strictly necessary since allocate will reset 
        # this array to empty before pushing it back on allocated stack. But this 
        # ensures the interval pool is in a consistent emptied state upon 
        # completion of this function.
        push!(a.unallocated, temp)
    end
end

# Allocate a zero length array which will be filled with Threads.nthreads() 
# entries by the __init__ method for OpticSim module. Have to do this in the 
#    __init__ method because this captures the load time environment. const 
#    values are evaluated at precompile time and the number of threads in these 
#    two environments can be different.
const threadedintervalpool = Vector{Dict{DataType, IntervalPool}}()
#const threadedintervalpool = [Dict{DataType,IntervalPool}(
# [Float64 => IntervalPool{Float64}()]) for _ in 1:Threads.nthreads()]

function newinintervalpool!(
    ::Type{T} = Float64,
    tid::Int = Threads.threadid(),
)::Vector{Interval{T}} where {T <: Real}
    if T ∉ keys(threadedintervalpool[tid])
        # if the type of the interval pool has changed then we need to refill 
        # it with the correct type
        threadedintervalpool[tid][T] = IntervalPool{T}()
    end
    return allocate!(threadedintervalpool[tid][T])
end

function indexednewinintervalpool!(
    ::Type{T} = Float64,
    tid::Int = Threads.threadid(),
)::Tuple{Int, Vector{Interval{T}}} where {T <: Real}
    a = newinintervalpool!(T, tid)
    index = length(threadedintervalpool[tid][T].allocated)
    return index, a
end

function emptyintervalpool!(
    ::Type{T} = Float64,
    tid::Int = Threads.threadid(),
) where {T}
    if T ∈ keys(threadedintervalpool[tid])
        empty!(threadedintervalpool[tid][T])
    end
end

getfromintervalpool(
    id::Int,
    tid::Int = Threads.threadid(),
)::Vector{Interval{Float64}} = getfromintervalpool(Float64, id, tid)

function getfromintervalpool(
    ::Type{T},
    id::Int,
    tid::Int = Threads.threadid(),
)::Vector{Interval{T}} where {T <: Real}
    return threadedintervalpool[tid][T].allocated[id]
end

# ============================================================================ #

macro inplaceinsertionsort(array, f)
    return esc(quote
        for i in 2:length($array)
            value = $array[i]
            j = i - 1
            while j > 0 && $f($array[j]) > $f(value)
                $array[j + 1] = $array[j]
                j = j - 1
            end
            $array[j + 1] = value
        end
    end)
end

"""
Datatype representing an ordered series of disjoint intervals on a ray.
An arbitrary array of `Interval`s can be input to the constructor and they will 
automatically be processed into a valid `DisjointUnion` 
(or a single [`Interval`](@ref) if appropriate).

```julia
DisjointUnion(intervals::AbstractVector{Interval{R}})
```
"""
struct DisjointUnion{T <: Real}
    intervalarrayindex::Int

    function DisjointUnion(
        a::Interval{R},
        b::Interval{R},
    )::Union{DisjointUnion{R}, Interval{R}} where {R <: Real}
        temp = newinintervalpool!(R)
        push!(temp, a)
        push!(temp, b)
        return DisjointUnion(temp)
    end

    function DisjointUnion(
        intervals::AbstractVector{Interval{R}},
    )::Union{DisjointUnion{R}, Interval{R}} where {R <: Real}
        if length(intervals) === 1
            return intervals[1]
        end

        @inplaceinsertionsort(intervals, upper)

        index, result = indexednewinintervalpool!(R)

        while length(intervals) > 0
            current = pop!(intervals)

            while length(intervals) > 0
                temp = intervalintersection(current, last(intervals))
                if !(temp isa EmptyInterval{R})
                    current = intervalunion(current, last(intervals))
                    pop!(intervals)
                else
                    break
                end
            end

            push!(result, current)
        end

        if length(result) === 1
            return result[1]
        else
            @inplaceinsertionsort(result, lower)
            return new{R}(index)
        end
    end
end

function intervals(a::DisjointUnion{R}) where {R <: Real}
    return getfromintervalpool(R, a.intervalarrayindex)
end
function Base.show(io::IO, a::DisjointUnion{R}) where {R <: Real}
    println(io, DisjointUnion{R})
    for i in intervals(a)
        show(io, i)
    end
end

Base.getindex(a::DisjointUnion, index::Int) = intervals(a)[index]
Base.firstindex(a::DisjointUnion) = firstindex(intervals(a))
Base.lastindex(a::DisjointUnion) = lastindex(intervals(a))
Base.length(a::DisjointUnion) = length(intervals(a))
function Base.:(==)(a::DisjointUnion, b::DisjointUnion)
    return intervals(a) == intervals(b)
end
# this works so long as the elements in the intervals array are structs 
# and not arrays.

# ============================================================================ #

"""
    isemptyinterval(a) -> Bool

Returns true if `a` is an [`EmptyInterval`](@ref). In performance critical 
contexts use `a isa EmptyInterval{T}`.

"""
isemptyinterval(::EmptyInterval) = true
isemptyinterval(::Interval) = false
isemptyinterval(::DisjointUnion) = false

"""
    ispositivehalfspace(a) -> Bool

Returns true if `upper(a)` is [`Infinity`](@ref). In performance critical 
contexts check directly i.e. `upper(a) isa Infinity{T}`.

"""
function ispositivehalfspace(a::Interval{T}) where {T <: Real}
    return upper(a) isa Infinity{T}
end

"""
    israyorigininterval(a) -> Bool

Returns true if `lower(a)` is [`RayOrigin`](@ref). In performance critical 
contexts check directly i.e. `lower(a) isa RayOrigin{T}`.

"""
function israyorigininterval(a::Interval{T}) where {T <: Real}
    return lower(a) isa RayOrigin{T}
end

function isinfiniteinterval(a::Interval{T}) where {T <: Real}
    return israyorigininterval(a) && ispositivehalfspace(a)
end

function positivehalfspace(a::RayOrigin{T}) where {T <: Real}
    return Interval(a, Infinity(T))
end
function positivehalfspace(a::Intersection{T, N}) where {T <: Real, N}
    return Interval(a, Infinity(T))
end

function rayorigininterval(a::Intersection{T, N}) where {T <: Real, N}
    return Interval(RayOrigin(T), a)
end # special kind of interval that has ray origin as lower bound. 
# α(a::RayOrigin) will always return 0 so CSG operations should work.
function rayorigininterval(a::Infinity{T}) where {T <: Real}
    return Interval(RayOrigin(T), a)
end

"""
    halfspaceintersection(a::Interval{T}) -> Intersection{T,3}

Returns the [`Intersection`](@ref) from a half space [`Interval`](@ref), 
throws an error if not a half space.
"""
function halfspaceintersection(
    a::Interval{T},
)::Intersection{T, 3} where {T <: Real}
    la = lower(a)
    ua = upper(a)
    if ua isa Infinity{T} && !(la isa RayOrigin{T})
        return la
    elseif la isa RayOrigin{T} && !(ua isa Infinity{T})
        return ua
    else
        return throw(ErrorException("Not a half-space: $a"))
    end
end

"""
    closestintersection(a::Union{EmptyInterval{T},Interval{T},DisjointUnion{T}},
    ignorenull::Bool = true) -> Union{Nothing,Intersection{T,3}}

Returns the closest [`Intersection`](@ref) from an [`Interval`](@ref) or 
[`DisjointUnion`](@ref). Ignores intersection with null interfaces if 
`ignorenull` is true. Will return `nothing` if there is no valid intersection.
"""
closestintersection(::EmptyInterval, ::Bool = true) = nothing
function closestintersection(
    a::Interval{T},
    ignorenull::Bool = true,
)::Union{Nothing, Intersection{T, 3}} where {T <: Real}
    la = lower(a)
    ua = upper(a)
    if la isa RayOrigin{T}
        if !(ua isa Infinity{T}) &&
           !(ignorenull && interface(ua) isa NullInterface{T})
            return ua
        else
            return nothing
        end
    elseif !(ignorenull && interface(la) isa NullInterface{T})
        return la
    else
        return nothing
    end
end
function closestintersection(
    a::DisjointUnion{T},
    ignorenull::Bool = true,
)::Union{Nothing, Intersection{T, 3}} where {T <: Real}
    for i in intervals(a)
        c = closestintersection(i, ignorenull)
        if c !== nothing
            return c
        end
    end
    return nothing
end

function reversenormal(a::Interval{T})::Interval{T} where {T <: Real}
    la = lower(a)
    ua = upper(a)
    if la isa RayOrigin{T}
        if ua isa Infinity{T}
            return a
        else
            return Interval(la, reversenormal(ua))
        end
    else
        if ua isa Infinity{T}
            return Interval(reversenormal(la), ua)
        else
            return Interval(reversenormal(la), reversenormal(ua))
        end
    end
end

function reversenormal(a::DisjointUnion{R})::DisjointUnion{R} where {R <: Real}
    intvls = newinintervalpool!(R)
    for i in intervals(a)
        push!(intvls, reversenormal(i))
    end
    return DisjointUnion(intvls)
end

# ============================================================================ #

macro intervalintersectionhigh(low)
    esc(quote
        if ua isa Infinity{R}
            if ub isa Infinity{R}
                return Interval($low, Infinity(R))
            else
                if ub <= $low
                    return EmptyInterval(R)
                end
                return Interval($low, ub)
            end
        else
            if ub isa Infinity{R}
                if ua <= $low
                    return EmptyInterval(R)
                end
                return Interval($low, ua)
            else
                high = min(ua, ub)
                if high <= $low
                    return EmptyInterval(R)
                end
                return Interval($low, high)
            end
        end
    end)
end

function intervalintersection(
    a::Interval{R},
    b::Interval{R},
)::Union{EmptyInterval{R}, Interval{R}} where {R <: Real}
    # This method is just doing the below, but to avoid type ambiguities 
    # things have to be much more complicated
    # if upper(a) <= lower(b) || upper(b) <= lower(a)
    #     return EmptyInterval(R)
    # end
    # low = max(lower(a), lower(b))
    # high = min(upper(a), upper(b))
    # return Interval(low, high)

    la = lower(a)
    lb = lower(b)
    ua = upper(a)
    ub = upper(b)
    if la isa RayOrigin{R}
        if lb isa RayOrigin{R}
            @intervalintersectionhigh(RayOrigin(R))
        else
            @intervalintersectionhigh(lb)
        end
    else
        if lb isa RayOrigin{R}
            @intervalintersectionhigh(la)
        else
            low = max(la, lb)
            @intervalintersectionhigh(low)
        end
    end
end

function intervalintersection(
    ::EmptyInterval{T},
    ::EmptyInterval{T},
) where {T <: Real}
    return EmptyInterval(T)
end
function intervalintersection(
    ::EmptyInterval{T},
    ::Interval{T},
) where {T <: Real}
    return EmptyInterval(T)
end
function intervalintersection(
    ::Interval{T},
    ::EmptyInterval{T},
) where {T <: Real}
    return EmptyInterval(T)
end
function intervalintersection(
    ::EmptyInterval{T},
    ::DisjointUnion{T},
) where {T <: Real}
    return EmptyInterval(T)
end
function intervalintersection(
    ::DisjointUnion{T},
    ::EmptyInterval{T},
) where {T <: Real}
    return EmptyInterval(T)
end
function intervalintersection(
    a::Interval{T},
    b::DisjointUnion{T},
) where {T <: Real}
    return intervalintersection(a, intervals(b))
end
function intervalintersection(
    a::DisjointUnion{T},
    b::Interval{T},
) where {T <: Real}
    return intervalintersection(b, intervals(a))
end
function intervalintersection(
    a::DisjointUnion{T},
    b::DisjointUnion{T},
) where {T <: Real}
    return intervalintersection(intervals(a), intervals(b))
end

function intervalintersection(
    a::Interval{T},
    b::AbstractVector{Interval{T}},
)::Union{EmptyInterval{T}, Interval{T}, DisjointUnion{T}} where {T <: Real}
    temp = nothing
    int1 = nothing
    for bint in b
        intsct = intervalintersection(a, bint)
        if !(intsct isa EmptyInterval{T})
            if int1 === nothing
                int1 = intsct
            elseif temp === nothing
                temp = newinintervalpool!(T)
                push!(temp, int1)
                push!(temp, intsct)
            else
                push!(temp, intsct)
            end
        end
    end
    if int1 === nothing
        return EmptyInterval(T)
    elseif int1 !== nothing && (temp === nothing)
        return int1
    else
        return DisjointUnion(temp)
    end
end

function intervalintersection(
    a::AbstractVector{Interval{T}},
    b::AbstractVector{Interval{T}},
)::Union{EmptyInterval{T}, Interval{T}, DisjointUnion{T}} where {T <: Real}
    temp = newinintervalpool!(T)
    for aint in a
        for bint in b
            intsct = intervalintersection(aint, bint)
            if !(intsct isa EmptyInterval{T})
                push!(temp, intsct)
            end
        end
    end
    if length(temp) == 0
        return EmptyInterval(T)
    else
        return DisjointUnion(temp)
    end
end

# ============================================================================ #

macro intervalunionhigh(low)
    esc(quote
        if ua isa Infinity{R}
            if ub isa Infinity{R}
                return Interval($low, Infinity(R))
            else
                if ub < la
                    return DisjointUnion(a, b)
                end
                return Interval($low, Infinity(R))
            end
        else
            if ub isa Infinity{R}
                if ua < lb
                    return DisjointUnion(a, b)
                end
                return Interval($low, Infinity(R))
            else
                if ua < lb || ub < la
                    return DisjointUnion(a, b)
                end
                return Interval($low, max(ua, ub))
            end
        end
    end)
end

function intervalunion(
    a::Interval{R},
    b::Interval{R},
)::Union{Interval{R}, DisjointUnion{R}} where {R <: Real}
    # This method is just doing the below, but to avoid type ambiguities things
    # have to be much more complicated
    # if upper(a) < lower(b) || upper(b) < lower(a)
    #     return DisjointUnion(a, b)
    # else
    #     low = min(lower(a), lower(b))
    #     high = max(upper(a), upper(b))
    #     return Interval(low, high)
    # end

    ua = upper(a)
    ub = upper(b)
    la = lower(a)
    lb = lower(b)
    if la isa RayOrigin{R}
        if lb isa RayOrigin{R}
            @intervalunionhigh(RayOrigin(R))
        else
            @intervalunionhigh(RayOrigin(R))
        end
    else
        if lb isa RayOrigin{R}
            @intervalunionhigh(RayOrigin(R))
        else
            low = min(la, lb)
            @intervalunionhigh(low)
        end
    end
end

function intervalunion(a::Interval{T}, b::DisjointUnion{T}) where {T <: Real}
    return intervalunion(b, a)
end
function intervalunion(a::DisjointUnion{T}, b::Interval{T}) where {T <: Real}
    temp = newinintervalpool!(T)
    for intvl in intervals(a)
        push!(temp, intvl)
    end
    push!(temp, b)
    return DisjointUnion(temp)
end
function intervalunion(
    a::DisjointUnion{T},
    b::DisjointUnion{T},
) where {T <: Real}
    temp = newinintervalpool!(T)
    for intvl in intervals(a)
        push!(temp, intvl)
    end
    for intvl in intervals(b)
        push!(temp, intvl)
    end
    return DisjointUnion(temp)
end
function intervalunion(::EmptyInterval{T}, ::EmptyInterval{T}) where {T <: Real}
    return EmptyInterval(T)
end
intervalunion(::EmptyInterval{T}, b::Interval{T}) where {T <: Real} = b
intervalunion(a::Interval{T}, ::EmptyInterval{T}) where {T <: Real} = a
intervalunion(::EmptyInterval{T}, b::DisjointUnion{T}) where {T <: Real} = b
intervalunion(a::DisjointUnion{T}, ::EmptyInterval{T}) where {T <: Real} = a

# ============================================================================ #

function intervalcomplement(a::DisjointUnion{T}) where {T <: Real}
    res = rayorigininterval(Infinity(T))
    for intrval in intervals(a)
        @assert !(intrval isa EmptyInterval{T}) //
                "should never have an empty interval in a DisjointUnion"
        compl =
            intervalcomplement(intrval)::Union{Interval{T}, DisjointUnion{T}}
        if compl isa Interval{T}
            res = intervalintersection(res, compl)
        else
            res = intervalintersection(res, compl)
        end
    end
    return res
end

function intervalcomplement(::EmptyInterval{T}) where {T <: Real}
    return Interval(RayOrigin(T), Infinity(T))
end
function intervalcomplement(a::Interval{T}) where {T <: Real}
    ua = upper(a)
    la = lower(a)
    if la isa RayOrigin{T}
        if ua isa Infinity{T}
            return EmptyInterval(T)
        else
            return positivehalfspace(ua) # this is not strictly correct
        end
    else
        if ua isa Infinity{T}
            return rayorigininterval(la)
        else
            return DisjointUnion(rayorigininterval(la), positivehalfspace(ua))
        end
    end
end

function difference(a::DisjointUnion, b::DisjointUnion)
    return intervalintersection(a, intervalcomplement(b))
end

"""
Apply a Transform to an Interval object
"""
function Base.:*(transformation::Transform{T}, a::Interval{T}) where {T <: Real}
    # looks ridiculous but necessary to dissambiguate the elements of the interval
    u = upper(a)
    l = lower(a)
    if l isa RayOrigin{T}
        if u isa Infinity{T}
            return Interval(l, u)
        else
            u = transformation * u
            return Interval(l, u)
        end
    else
        l = transformation * l
        if u isa Infinity{T}
            return Interval(l, u)

        else
            u = transformation * u
            return Interval(l, u)
        end
    end
end

end # module
