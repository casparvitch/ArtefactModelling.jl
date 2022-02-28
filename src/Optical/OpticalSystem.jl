# Copied/adapted from code under:
# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    AbstractOpticalSystem{T<:Real}

Abstract type for any optical system, must parameterized by the datatype of 
entities within the system `T` (numeric).

The OpticalSystem will hold relevant global parameters, as well as holding
all of the objects in the scene.
"""
abstract type AbstractOpticalSystem{T <: Real} end

# ============================================================================ #
# ====                        CSGOpticalSystem                            ==== #
# ============================================================================ #

"""
    CSGOpticalSystem{T,D<:Real,S<:Surface{T},L<:LensAssembly{T}} <: AbstractOpticalSystem{T}

An optical system containing a lens assembly with all optical elements and a 
detector surface with associated image. The system can be at a specified 
temperature and pressure.

There are two number types in the type signature. 
The `T` type parameter is the numeric type for geometry in the optical system, 
the `D` type parameter is the numeric type of the pixels in the detector image. 
This way you can have `Float64` geometry, where high precision is essential, 
but the pixels in the detector can be `Float32` since precision is much less 
critical for image data, or Complex if doing wave optic simulations.

The detector can be any [`Surface`](@ref) 
which implements [`uv`](@ref), [`uvtopix`](@ref) and [`onsurface`](@ref),
typically this is one of 
[`Rectangle`](@ref), [`Ellipse`](@ref) or [`SphericalCap`](@ref).

```julia
CSGOpticalSystem(
    assembly::LensAssembly,
    detector::Surface,
    detectorpixelsx = 1000,
    detectorpixelsy = 1000, ::Type{D} = Float32;
    temperature = OpticSim.GlassCat.TEMP_REF,
    pressure = OpticSim.GlassCat.PRESSURE_REF
)
```
"""
struct CSGOpticalSystem{
    T,
    D <: Number,
    S <: Surface{T}, # should allow for more than one detector?
    L <: LensAssembly{T}, # need a new type for 'objects_in_scene' type thing
} <: AbstractOpticalSystem{T}
    assembly::L
    detector::S
    detectorimage::HierarchicalImage{D}
    temperature::T
    pressure::T

    function CSGOpticalSystem(
        assembly::L,
        detector::S,
        detectorpixelsx::Int = 1000,
        detectorpixelsy::Int = 1000,
        ::Type{D} = Float32;
        temperature::Union{T, Unitful.Temperature} = convert(T, TEMP_REF),
        pressure::T = convert(T, PRESSURE_REF),
    ) where {T <: Real, S <: Surface{T}, L <: LensAssembly{T}, D <: Number}
        @assert hasmethod(uv, (S, SVector{3, T})) "Detector must implement uv()"
        @assert hasmethod(uvtopix, (S, SVector{2, T}, Tuple{Int, Int})) "Detector must implement uvtopix()"
        @assert hasmethod(onsurface, (S, SVector{3, T})) "Detector must implement onsurface()"
        opticalinterface = interface(detector)
        @assert insidematerialid(opticalinterface) ==
                outsidematerialid(opticalinterface) "Detector must have same material either side"
        @assert interface(detector) !== NullInterface(T) "Detector can't have null interface"
        image = HierarchicalImage{D}(detectorpixelsy, detectorpixelsx)
        if temperature isa Unitful.Temperature
            temperature = Unitful.ustrip(T, °C, temperature)
        end
        return new{T, D, S, L}(
            assembly,
            detector,
            image,
            temperature,
            convert(T, pressure),
        )
    end
end

function Base.copy(a::CSGOpticalSystem)
    return CSGOpticalSystem(
        a.assembly,
        a.detector,
        size(a.detectorimage)...,
        temperature = (a.temperature)°C,
        pressure = a.pressure,
    )
end

# added this show method because the type of CSGOpticalSystem is gigantic and printing it in the REPL can crash the
# system
function Base.show(io::IO, a::CSGOpticalSystem{T}) where {T}
    return print(
        io,
        "CSGOpticalSystem{$T}($(temperature(a))," *
        " $(pressure(a)), $(assembly(a)), $(detector(a)))",
    )
end

"""
    assembly(system::AbstractOpticalSystem{T}) -> LensAssembly{T}

Get the [`LensAssembly`](@ref) of `system`.
"""
assembly(system::CSGOpticalSystem{T}) where {T <: Real} = system.assembly

detector(system::CSGOpticalSystem) = system.detector

"""
    detectorimage(system::AbstractOpticalSystem{T}) -> HierarchicalImage{D}

Get the detector image of `system`.
`D` is the datatype of the detector image and is not necessarily the same as 
the datatype of the system `T`.
"""
detectorimage(system::CSGOpticalSystem) = system.detectorimage

detectorsize(system::CSGOpticalSystem) = size(system.detectorimage)

"""
    temperature(system::AbstractOpticalSystem{T}) -> T

Get the temperature of `system` in °C.
"""
temperature(system::CSGOpticalSystem{T}) where {T <: Real} = system.temperature

"""
    pressure(system::AbstractOpticalSystem{T}) -> T

Get the pressure of `system` in Atm.
"""
pressure(system::CSGOpticalSystem{T}) where {T <: Real} = system.pressure
"""
    resetdetector!(system::AbstractOpticalSystem{T})

Reset the deterctor image of `system` to zero.
"""
function resetdetector!(system::CSGOpticalSystem{T}) where {T <: Real}
    return reset!(system.detectorimage)
end

# not sure what the below does
# Base.Float32(a::T) where {T<:ForwardDiff.Dual} = Float32(ForwardDiff.value(a))

# ============================================================================ #
# ====                              Trace                                 ==== #
# ============================================================================ #

"""
    trace(system::AbstractOpticalSystem{T}, ray::OpticalRay{T}; trackrays = nothing, test = false)

Traces `system` with `ray`, if `test` is enabled then fresnel reflections are disabled and the power distribution will
not be correct. Returns either a [`LensTrace`](@ref) if the ray hits the detector or `nothing` otherwise.

`trackrays` can be passed an empty vector to accumulate the `LensTrace` objects at each intersection of `ray` with a
surface in the system.
"""
function trace(
    system::CSGOpticalSystem{T, D},
    r::OpticalRay{T, N};
    trackrays::Union{Nothing, Vector{LensTrace{T, N}}} = nothing,
    test::Bool = false,
) where {T <: Real, N, D <: Number}
    if power(r) < POWER_THRESHOLD
        return nothing
    end

    result = trace(
        system.assembly,
        r,
        temperature(system),
        pressure(system),
        trackrays = trackrays,
        test = test,
    )

    if result === nothing || result === nopower
        emptyintervalpool!(T)
        return nothing
    else #ray intersected lens assembly so continue to see if ray intersects detector
        intsct = surfaceintersection(detector(system), ray(result))
        if intsct === nothing # no intersection of final ray with detector
            emptyintervalpool!(T)
            return nothing
        end

        detintsct = closestintersection(intsct)
        if detintsct === nothing
            emptyintervalpool!(T)
            return nothing
        else
            # need to modify power and path length accordingly for the intersection with the detector
            surfintsct = point(detintsct)
            nml = normal(detintsct)
            opticalinterface = interface(detintsct)::FresnelInterface{T}
            λ = wavelength(r)

            # Optical path length is measured in mm
            # in this case the result ray is exact so no correction for RAY_OFFSET is needed
            geometricpathlength = norm(surfintsct - origin(ray(result)))
            opticalpathlength = geometricpathlength
            pow = power(result)

            m = outsidematerialid(opticalinterface)
            # compute updated power based on absorption coefficient of material using Beer's law
            # this will almost always not apply as the detector will be in air, but it's possible that the detector is
            # not in air, in which case this is necessary
            if !isair(m)
                mat::Glass = glassforid(m)
                nᵢ = index(
                    mat,
                    λ,
                    temperature = temperature(system),
                    pressure = pressure(system),
                )::T
                α = absorption(
                    mat,
                    λ,
                    temperature = temperature(system),
                    pressure = pressure(system),
                )::T
                if α > zero(T)
                    internal_trans = exp(-α * geometricpathlength)
                    if rand() >= internal_trans
                        return nothing
                    end
                    pow = pow * internal_trans
                end
                opticalpathlength = nᵢ * geometricpathlength
            end

            temp = LensTrace{T, N}(
                OpticalRay(
                    ray(ray(result)),
                    pow,
                    wavelength(result),
                    opl = pathlength(result) + opticalpathlength,
                    nhits = nhits(result) + 1,
                    sourcenum = sourcenum(r),
                    sourcepower = sourcepower(r),
                ),
                detintsct,
            )
            if trackrays !== nothing
                push!(trackrays, temp)
            end

            # increment the detector image
            pixu, pixv = uvtopix(
                detector(system),
                uv(detintsct),
                size(system.detectorimage),
            )
            system.detectorimage[pixv, pixu] += convert(D, sourcepower(r)) # TODO will need to handle different detector
            #      image types a bit better than this

            # should be okay to assume intersection will not be a DisjointUnion for all the types of detectors we will
            # be using
            emptyintervalpool!(T)
            return temp
        end
    end
end
