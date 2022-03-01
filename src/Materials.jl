
module Materials

using StaticArrays: SVector

abstract type AbstractMaterial end

abstract type AbstractTransparentMaterial <: AbstractMaterial end
# add absorption etc. methods, they'll be useful

# ============================================================================ #
# ====                               Air                                  ==== #
# ============================================================================ #

struct Air <: AbstractTransparentMaterial end

name(::Air) = :Air

function index(
    ::Air,
    ::T;
    temperature::T = T(20.0),
    pressure::T = T(1.0),
)::T where {T <: Real}
    return one(T)
end

function absorption(
    ::Air,
    ::T;
    temperature::T = T(20.0),
    pressure::T = T(1.0),
)::T where {T <: Real}
    return zero(T)
end

"""λ in m."""
function absairindex(
    λ::T;
    temperature::T = T(20.0),
    pressure::T = T(1.0),
)::T where {T <: Real}
    # convert to required units
    λ *= 1e-6 # convert to um
    n_ref =
        one(T) + (
            (
                6432.8 +
                ((2949810.0 * λ^2) / (146.0 * λ^2 - one(T))) +
                ((25540.0 * λ^2) / (41.0 * λ^2 - one(T)))
            ) * 1e-8
        )
    n_rel =
        one(T) +
        ((n_ref - one(T)) / (one(T) + (temperature - 15.0) * 0.0034785)) *
        pressure
    return n_rel
end

# assume wavelengths are in m
function polyfit_indices(
    wavelengths::Union{AbstractRange{<:Real}, AbstractArray{<:Real, 1}},
    indices::AbstractArray{<:Number, 1};
    degree::Int = 5,
)
    w = 1e-6 .* wavelengths # convert to um
    okay = (indices .> 0.0)
    if !any(okay)
        return (ones(Float64, size(w)) .* NaN, nothing)
    end
    xs = range(-1.0, stop = 1.0, length = length(w[okay]))
    poly = fit(xs, indices[okay], degree)
    interp_indices = poly.(xs)
    # ensure output has all entries
    out = ones(Float64, size(w)) .* NaN
    out[okay] = interp_indices
    return (out, poly)
end

# ============================================================================ #
# ====                          Perfect Abs.                              ==== #
# ============================================================================ #

struct PerfectAbsorber <: AbstractMaterial end

name(::PerfectAbsorber) = :PerfectAbsorber

function absorption(
    ::M,
    ::T;
    temperature::T = T(20.0),
    pressure::T = T(1.0),
)::T where {M <: AbstractMaterial, T <: Real}
    return one(T)
end

# ============================================================================ #
# ====                             Glass                                  ==== #
# ============================================================================ #

struct Glass <: AbstractTransparentMaterial
    name::Symbol
    dispform::Int
    C1::Float64
    C2::Float64
    C3::Float64
    C4::Float64
    C5::Float64
    C6::Float64
    C7::Float64
    C8::Float64
    C9::Float64
    C10::Float64
    λmin::Float64
    λmax::Float64
    D₀::Float64
    D₁::Float64
    D₂::Float64
    E₀::Float64
    E₁::Float64
    λₜₖ::Float64
    temp::Float64
    ΔPgF::Float64
    PR::Float64
    relcost::Float64
    TCE::Float64
    CR::Float64
    status::Int
    SR::Float64
    transmission::Union{Nothing, SVector{100, SVector{3, Float64}}}
    transmissionN::Int
    Nd::Float64
    AR::Float64
    FR::Float64
    exclude_sub::Int
    Vd::Float64
    ignore_thermal_exp::Int
    p::Float64
    meltfreq::Int
    function Glass(
        name::Symbol,
        dispform,
        C1,
        C2,
        C3,
        C4,
        C5,
        C6,
        C7,
        C8,
        C9,
        C10,
        λmin,
        λmax,
        D₀,
        D₁,
        D₂,
        E₀,
        E₁,
        λₜₖ,
        temp,
        ΔPgF,
        PR,
        relcost,
        TCE,
        CR,
        status,
        SR,
        transmission,
        Nd,
        AR,
        FR,
        exclude_sub,
        Vd,
        ignore_thermal_exp,
        p,
        meltfreq,
    )
        # need a constructor to massage the transmission data
        if transmission === nothing
            transmission_s = nothing
            transmissionN = -1
        else
            fill =
                [SVector(0.0, 0.0, 0.0) for _ in 1:(100 - length(transmission))]
            transmission_s =
                SVector{100, SVector{3, Float64}}(transmission..., fill...)
            transmissionN = length(transmission)
        end
        return new(
            name,
            dispform,
            C1,
            C2,
            C3,
            C4,
            C5,
            C6,
            C7,
            C8,
            C9,
            C10,
            λmin,
            λmax,
            D₀,
            D₁,
            D₂,
            E₀,
            E₁,
            λₜₖ,
            temp,
            ΔPgF,
            PR,
            relcost,
            TCE,
            CR,
            status,
            SR,
            transmission_s,
            transmissionN,
            Nd,
            AR,
            FR,
            exclude_sub,
            Vd,
            ignore_thermal_exp,
            p,
            meltfreq,
        )
    end
end

name(g::Glass) = g.name

"""
    info([io::IO], glass::AbstractGlass)

Print out all data associated with `glass` in an easily readable format.

# Examples
```julia-repl
julia> info(GlassCat.RPO.IG4)
name:                                              :AGF_52
Dispersion formula:                                Schott (1)
Dispersion formula coefficients:
     a₀:                                           6.91189161
     a₁:                                           -0.000787956404
     a₂:                                           -4.22296071
     a₃:                                           142.900646
     a₄:                                           -1812.32748
     a₅:                                           7766.33028
Valid wavelengths:                                 3.0μm to 12.0μm
Reference temperature:                              20.0°C
Thermal ΔRI coefficients:
     D₀:                                           3.24e-5
     D₁:                                           0.0
     D₂:                                           0.0
     E₀:                                           0.0
     E₁:                                           0.0
     λₜₖ:                                          0.0
TCE (÷1e-6):                                       20.4
Ignore thermal expansion:                          false
Density (p):                                       4.47g/m³
ΔPgF:                                              0.0
RI at sodium D-Line (587nm):                       1.0
Abbe Number:                                       0.0
Cost relative to N_BK7:                              ?
Status:                                            Standard (0)
Melt frequency:                                    0
Exclude substitution:                              false
```
"""
function info(io::IO, glass::Glass)
    D = glass.dispform

    println(io, "$(rpad("name:", 48)) $(glass.name)")

    if (D == -2)
        println(io, "$(rpad("Dispersion formula:", 50)) Cauchy (-2)")
    elseif (D == -1)
        println(
            io,
            "$(rpad("Dispersion formula:", 50)) Fitted for model/MIL glass",
        )
    else
        println(
            io,
            "$(rpad("Dispersion formula:", 50)) $(DISPFORM_NAMES[D]) ($D)",
        )
    end
    println(io, "Dispersion formula coefficients:")
    if (D == -2) # Cauchy
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C:", 45)) $(glass.C3)")
        println(io, "     $(rpad("D:", 45)) $(glass.C4)")
        println(io, "     $(rpad("E:", 45)) $(glass.C5)")
        println(io, "     $(rpad("F:", 45)) $(glass.C6)")
    elseif (D == -1)
        println(io, "     $(rpad("C₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("C₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("C₃:", 45)) $(glass.C4)")
    elseif (D == 1) # Schott
        println(io, "     $(rpad("a₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("a₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("a₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("a₃:", 45)) $(glass.C4)")
        println(io, "     $(rpad("a₄:", 45)) $(glass.C5)")
        println(io, "     $(rpad("a₅:", 45)) $(glass.C6)")
    elseif (D == 2)  # Sellmeier1
        println(io, "     $(rpad("K₁:", 45)) $(glass.C1)")
        println(io, "     $(rpad("L₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("K₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("L₂:", 45)) $(glass.C4)")
        println(io, "     $(rpad("K₃:", 45)) $(glass.C5)")
        println(io, "     $(rpad("L₃:", 45)) $(glass.C6)")
    elseif (D == 3)  # Herzberger
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C:", 45)) $(glass.C3)")
        println(io, "     $(rpad("D:", 45)) $(glass.C4)")
        println(io, "     $(rpad("E:", 45)) $(glass.C5)")
        println(io, "     $(rpad("F:", 45)) $(glass.C6)")
    elseif (D == 4)  # Sellmeier2
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("λ₁:", 45)) $(glass.C3)")
        println(io, "     $(rpad("B₂:", 45)) $(glass.C4)")
        println(io, "     $(rpad("λ₂:", 45)) $(glass.C5)")
    elseif (D == 5)  # Conrady
        println(io, "     $(rpad("n₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("A:", 45)) $(glass.C2)")
        println(io, "     $(rpad("B:", 45)) $(glass.C3)")
    elseif (D == 6)  # Sellmeier3
        println(io, "     $(rpad("K₁:", 45)) $(glass.C1)")
        println(io, "     $(rpad("L₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("K₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("L₂:", 45)) $(glass.C4)")
        println(io, "     $(rpad("K₃:", 45)) $(glass.C5)")
        println(io, "     $(rpad("L₃:", 45)) $(glass.C6)")
        println(io, "     $(rpad("K₄:", 45)) $(glass.C7)")
        println(io, "     $(rpad("L₄:", 45)) $(glass.C8)")
    elseif (D == 7) || (D == 8)  # HandbookOfOptics1/2
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C:", 45)) $(glass.C3)")
        println(io, "     $(rpad("D:", 45)) $(glass.C4)")
    elseif (D == 9)  # Sellmeier4
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C:", 45)) $(glass.C3)")
        println(io, "     $(rpad("D:", 45)) $(glass.C4)")
        println(io, "     $(rpad("E:", 45)) $(glass.C5)")
    elseif (D == 10) || (D == 12) # Extended1/2
        println(io, "     $(rpad("a₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("a₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("a₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("a₃:", 45)) $(glass.C4)")
        println(io, "     $(rpad("a₄:", 45)) $(glass.C5)")
        println(io, "     $(rpad("a₅:", 45)) $(glass.C6)")
        println(io, "     $(rpad("a₆:", 45)) $(glass.C7)")
        println(io, "     $(rpad("a₇:", 45)) $(glass.C8)")
    elseif (D == 11)  # Sellmeier5
        println(io, "     $(rpad("K₁:", 45)) $(glass.C1)")
        println(io, "     $(rpad("L₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("K₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("L₂:", 45)) $(glass.C4)")
        println(io, "     $(rpad("K₃:", 45)) $(glass.C5)")
        println(io, "     $(rpad("L₃:", 45)) $(glass.C6)")
        println(io, "     $(rpad("K₄:", 45)) $(glass.C7)")
        println(io, "     $(rpad("L₄:", 45)) $(glass.C8)")
        println(io, "     $(rpad("K₅:", 45)) $(glass.C9)")
        println(io, "     $(rpad("L₅:", 45)) $(glass.C10)")
    elseif (D == 13)  # Extended3
        println(io, "     $(rpad("a₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("a₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("a₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("a₃:", 45)) $(glass.C4)")
        println(io, "     $(rpad("a₄:", 45)) $(glass.C5)")
        println(io, "     $(rpad("a₅:", 45)) $(glass.C6)")
        println(io, "     $(rpad("a₆:", 45)) $(glass.C7)")
        println(io, "     $(rpad("a₇:", 45)) $(glass.C8)")
        println(io, "     $(rpad("a₈:", 45)) $(glass.C9)")
    else
        println(io, "     INVALID DISPERSION FORMULA!!")
    end
    println(
        io,
        "$(rpad("Valid wavelengths:", 50)) $(glass.λmin)μm to $(glass.λmax)μm",
    )
    println(io, "$(rpad("Reference temperature:", 50)) $(glass.temp)°C")

    if !isnan(glass.D₀) && (
        glass.D₀ != 0 ||
        glass.D₁ != 0 ||
        glass.D₂ != 0 ||
        glass.E₀ != 0 ||
        glass.E₁ != 0
    )
        println(io, "Thermal ΔRI coefficients:")
        println(io, "     $(rpad("D₀:", 45)) $(glass.D₀)")
        println(io, "     $(rpad("D₁:", 45)) $(glass.D₁)")
        println(io, "     $(rpad("D₂:", 45)) $(glass.D₂)")
        println(io, "     $(rpad("E₀:", 45)) $(glass.E₀)")
        println(io, "     $(rpad("E₁:", 45)) $(glass.E₁)")
        println(io, "     $(rpad("λₜₖ:", 45)) $(glass.λₜₖ)")
    end

    println(io, "$(rpad("TCE (÷1e-6):", 50)) $(glass.TCE)")
    println(
        io,
        "$(rpad("Ignore thermal expansion:", 50)) $(glass.ignore_thermal_exp == 1)",
    )

    println(io, "$(rpad("Density (p):", 50)) $(glass.p)g/m³")
    println(io, "$(rpad("ΔPgF:", 50)) $(glass.ΔPgF)")

    println(io, "$(rpad("RI at sodium D-Line (587nm):", 50)) $(glass.Nd)")
    println(io, "$(rpad("Abbe Number:", 50)) $(glass.Vd)")

    println(
        io,
        "$(rpad("Cost relative to N_BK7:", 50)) $(glass.relcost == -1 ? "?" : glass.relcost)",
    )

    if glass.CR != -1 ||
       glass.FR != -1 ||
       glass.SR != -1 ||
       glass.AR != -1 ||
       glass.PR != -1
        println(io, "Environmental resistance:")
        println(
            io,
            "     $(rpad("Climate (CR):", 45)) $(glass.CR == -1 ? "?" : glass.CR)",
        )
        println(
            io,
            "     $(rpad("Stain (FR):", 45)) $(glass.FR == -1 ? "?" : glass.FR)",
        )
        println(
            io,
            "     $(rpad("Acid (SR):", 45)) $(glass.SR == -1 ? "?" : glass.SR)",
        )
        println(
            io,
            "     $(rpad("Alkaline (AR):", 45)) $(glass.AR == -1 ? "?" : glass.AR)",
        )
        println(
            io,
            "     $(rpad("Phosphate (PR):", 45)) $(glass.PR == -1 ? "?" : glass.PR)",
        )
    end

    println(
        io,
        "$(rpad("Status:", 50)) $(STATUS[glass.status + 1]) ($(glass.status))",
    )
    println(
        io,
        "$(rpad("Melt frequency:", 50)) $(glass.meltfreq == -1 ? "?" : glass.meltfreq)",
    )
    println(
        io,
        "$(rpad("Exclude substitution:", 50)) $(glass.exclude_sub == 1)",
    )

    if glass.transmission !== nothing
        println(io, "Transmission data:")
        println(
            io,
            "$(lpad("Wavelength", 15))$(lpad("Transmission", 15))$(lpad("Thickness", 15))",
        )
        for i in 1:(glass.transmissionN)
            λ, t, τ = glass.transmission[i]
            println(
                io,
                "$(lpad("$(λ)μm", 15))$(lpad(t, 15))$(lpad("$(τ)mm", 15))",
            )
        end
    end
end

NBK7 = Glass(
    :NBK7,
    2,
    1.03961212,
    0.00600069867,
    0.231792344,
    0.0200179144,
    1.01046945,
    103.560653,
    0.0,
    0.0,
    NaN,
    NaN,
    0.3,
    2.5,
    1.86e-6,
    1.31e-8,
    -1.37e-11,
    4.34e-7,
    6.27e-10,
    0.17,
    20.0,
    -0.0009,
    2.3,
    1.0,
    7.1,
    1.0,
    1,
    1.0,
    [
        (0.3, 0.05, 25.0),
        (0.31, 0.25, 25.0),
        (0.32, 0.52, 25.0),
        (0.334, 0.78, 25.0),
        (0.35, 0.92, 25.0),
        (0.365, 0.971, 25.0),
        (0.37, 0.977, 25.0),
        (0.38, 0.983, 25.0),
        (0.39, 0.989, 25.0),
        (0.4, 0.992, 25.0),
        (0.405, 0.993, 25.0),
        (0.42, 0.993, 25.0),
        (0.436, 0.992, 25.0),
        (0.46, 0.993, 25.0),
        (0.5, 0.994, 25.0),
        (0.546, 0.996, 25.0),
        (0.58, 0.995, 25.0),
        (0.62, 0.994, 25.0),
        (0.66, 0.994, 25.0),
        (0.7, 0.996, 25.0),
        (1.06, 0.997, 25.0),
        (1.53, 0.98, 25.0),
        (1.97, 0.84, 25.0),
        (2.325, 0.56, 25.0),
        (2.5, 0.36, 25.0),
    ],
    1.5168,
    2.3,
    0.0,
    0,
    64.17,
    0,
    2.51,
    0.0,
)

"""
    modelglass(Nd::Float64, Vd::Float64, ΔPgF::Float64, name::Symbol) -> Glass

Generates a glass object for the given refractive index at d-light (587.5618nm),
 `Nd`, the Abbe number also at d-light, `Vd`, and partial dispersion, `ΔPgF`.
The mean error to measured data for these models is typically small - 
usually < 0.0001.
Behavior may differ from other optical simulation tools when using model glasses.

The approximate dispersion calculation used for these glasses is generally only 
valid for visible wavelengths, in this case a limit of 360nm to 750nm is imposed.

# Examples
```julia-repl
julia> index(modelglass(1.5168, 64.17, 0.0), 0.5875618)
1.5168003970108495

julia> index(modelglass(1.2344, 61.57, 0.003), 0.678)
1.2329425902693352
```
"""
function modelglass(Nd::Float64, Vd::Float64, ΔPgF::Float64, name::Symbol)
    # from Schott "TIE-29: Refractive Index and Dispersion"
    a = ΔPgF + 0.6438 - 0.001682 * Vd
    # Using fitting results from https://www.gnu.org/software/goptical/manual/Material_Abbe_class_reference.html
    C1 = a * -6.11873891971188577088 + 1.17752614766485175224
    C2 = a * 18.27315722388047447566 + -8.93204522498095698779
    C3 = a * -14.55275321129051135927 + 7.91015964461522003148
    C4 = a * 3.48385106908642905310 + -1.80321117937358499361
    return Glass(
        name,
        -1,
        C1,
        C2,
        C3,
        C4,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.36,
        0.75,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        TEMP_REF,
        ΔPgF,
        -1.0,
        -1.0,
        0.0,
        -1.0,
        0,
        -1.0,
        nothing,
        Nd,
        -1.0,
        -1.0,
        0,
        Vd,
        1,
        0.0,
        0.0,
    )
end

"""
    glassfromMIL(glasscode::Union{Float64,Int}, name::Symbol) -> Glass

Generates a glass object for the given glass code based on U.S. military 
standard MIL-G-174, see 
[the MIL specification](http://www.newportglass.com/GeneCd.htm) 
for further details.

The glass code is a six-digit number specifying the glass according to its 
refractive index `Nd` at d-light (587.5618nm), and its Abbe number `Vd` also 
taken at d-light.
The resulting glass code is the value of `Nd - 1` rounded to three digits, 
followed by `Vd` rounded to three digits, with all decimal points ignored.
For example, N_BK7 has `Nd = 1.5168` and `Vd = 64.17`, giving a six-digit glass
code of `517642`.

For `Nd > 1.999` the format `1.123642` can be used representing `Nd = 2.123` and 
`Vd = 64.2`.

**Accuracy is poor given the low precision of the input parameters**, the mean 
error to measured data may be significant.
Behavior may differ from other optical simulation tools when using MIL glasses.
The approximate dispersion calculation used these glasses is generally only 
valid  for visible wavelengths, in this case a limit of 360nm to 750nm is 
imposed.

# Examples
```julia-repl
julia> index(glassfromMIL(517642), 0.5875618)
1.5170003960064509

julia> index(glassfromMIL(1.134642), 0.5875618)
2.1340008686098946
```
"""
function glassfromMIL(glasscode::Int, name::Symbol)
    Nd = floor(Int, glasscode / 1000) / 1000 + 1
    Vd = (glasscode - floor(Int, glasscode / 1000) * 1000) / 10
    g = modelglass(Nd, Vd, 0.0, name)
    return g
end

function glassfromMIL(glasscode::Float64, name::Symbol)
    @assert glasscode > 1.0
    glasscodeid = round(Int, glasscode * 1000000)
    Nd = floor(Int, glasscode * 1000) / 1000 + 1
    Vd = round(
        (glasscode * 1000 - floor(Int, glasscode * 1000)) * 100,
        digits = 1,
    )
    g = modelglass(Nd, Vd, 0.0, name)
    return g
end

info(g::Glass) = info(stdout, g)

"""
    absorption(glass::Glass, wavelength; temperature=20.0, pressure=1.0)

Compute the intensity absorption per mm of `glass` at `wavelength`, 
optionally at specified `temperature` and `pressure`.
Transmission values are linearly interpolated from the adjacent values in the 
data table of `glass`, if `wavelength` is below the minimum or above the 
maximum in the table then the nearest value is taken.

Absorption is defined as ``\\frac{-\\log(t)}{\\tau}`` where ``t`` is the 
transmission value and ``\\tau`` is the thickness, both of which are provided 
in the data table.

Wavelength, temp and pressure in m, °C and Atm. 
"""
function absorption(
    glass::Glass,
    λ::T;
    temperature::T = T(20.0),
    pressure::T = T(1.0),
)::T where {T <: Real}
    # if the glass has no transmission data then assume no absorption
    if glass.transmission === nothing
        return zero(T)
    end

    reference_temp = T(glass.temp)

    # to work out the wavelength at the reference temperature we need the RIs of air at system temp and at reference temp
    n_air_at_sys =
        absairindex(λ, temperature = temperature, pressure = pressure)
    n_air_at_ref = absairindex(λ, temperature = reference_temp)

    # scale the wavelength to air at the reference temperature/pressure (and um)
    λ = λ * 1e-6 * (n_air_at_sys / n_air_at_ref)

    tdata = glass.transmission
    N = glass.transmissionN
    if λ < tdata[1][1]
        t = tdata[1][2]
        τ = tdata[1][3]
        return T(-log1p(t - 1.0) / τ)
    elseif λ > tdata[N][1]
        t = tdata[N][2]
        τ = tdata[N][3]
        return T(-log1p(t - 1.0) / τ)
    else
        let λlow = 0.0,
            tlow = 0.0,
            τlow = 0.0,
            λhigh = 0.0,
            thigh = 0.0,
            τhigh = 0.0

            for i in 2:N
                if λ <= tdata[i][1]
                    λlow, tlow, τlow = tdata[i - 1]
                    λhigh, thigh, τhigh = tdata[i]
                    break
                end
            end
            λhigh = T(λhigh)
            λlow = T(λlow)
            δλ = λhigh - λlow
            @assert τlow == τhigh
            t = (tlow * (λhigh - λ) / δλ) + (thigh * (λ - λlow) / δλ)
            return -log1p(t - 1.0) / τhigh
        end
    end
end

"""
    index(glass::AbstractGlass, wavelength; temperature=20°C, pressure=1Atm)

Compute the refractive index of `glass` at `wavelength`, optionally at specified `temperature` and `pressure`.
Result is relative to the refractive index of air at given temperature and pressure.

Arguments are interpretted as m, °C and Atm respectively.

**This is defined to always equal 1.0 for Air at any temperature and pressure**,
 use [`absairindex`](@ref) for the absolute refractive index of air at a given 
 temperature and pressure.
"""
function index(
    glass::Glass,
    λ::T;
    temperature::T = T(20.0),
    pressure::T = T(1.0),
)::T where {T <: Real}
    # all calculations for the material must be done at the refernce temperature
    reference_temp = T(glass.temp)

    # to work out the wavelength at the reference temperature we need the RIs of air at system temp and at reference temp
    n_air_at_sys =
        absairindex(λ, temperature = temperature, pressure = pressure)
    n_air_at_ref = absairindex(λ, temperature = reference_temp)

    # scale the wavelength to air at the reference temperature/pressure
    λabs = λ * n_air_at_sys * 1e-6 # (and um)
    λ = λabs / n_air_at_ref

    if (λ < glass.λmin) || (λ > glass.λmax)
        error(
            "Cannot calculate an index for the specified wavelength: $λ, valid range: [$(glass.λmin), $(glass.λmax)].\n",
        )
    end

    if glass.dispform == -2
        # Cauchy
        n_rel =
            T(glass.C1) +
            (glass.C2 * λ^(-2)) +
            (glass.C3 * λ^(-4)) +
            (glass.C4 * λ^(-6)) +
            (glass.C5 * λ^(-8)) +
            (glass.C6 * λ^(-10))
    elseif glass.dispform == -1
        # use fitted result from GOptical:
        n_rel = T(
            glass.Nd +
            (glass.Nd - one(T)) / glass.Vd *
            (glass.C1 + glass.C2 / λ + glass.C3 / λ^2 + glass.C4 / λ^3),
        )
    elseif glass.dispform == 1
        # Schott
        formula_rhs =
            T(glass.C1) +
            (glass.C2 * λ^2) +
            (glass.C3 * λ^(-2)) +
            (glass.C4 * λ^(-4)) +
            (glass.C5 * λ^(-6)) +
            (glass.C6 * λ^(-8))
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 2
        # Sellmeier1
        formula_rhs =
            (glass.C1 * λ^2 / (λ^2 - glass.C2)) +
            (glass.C3 * λ^2 / (λ^2 - glass.C4)) +
            (glass.C5 * λ^2 / (λ^2 - glass.C6))
        n_rel = sqrt(formula_rhs + one(T))
    elseif glass.dispform == 3
        # Herzberger
        L = one(T) / (λ^2 - T(0.028))
        n_rel =
            T(glass.C1) +
            (glass.C2 * L) +
            (glass.C3 * L^2) +
            (glass.C4 * λ^2) +
            (glass.C5 * λ^4) +
            (glass.C6 * λ^6)
    elseif glass.dispform == 4
        # Sellmeier2
        formula_rhs =
            T(glass.C1) +
            (glass.C2 * λ^2 / (λ^2 - (glass.C3)^2)) +
            (glass.C4 * λ^2 / (λ^2 - (glass.C5)^2))
        n_rel = sqrt(formula_rhs + one(T))
    elseif glass.dispform == 5
        # Conrady
        n_rel = T(glass.C1) + (glass.C2 / λ) + (glass.C3 / λ^3.5)
    elseif glass.dispform == 6
        # Sellmeier3
        formula_rhs =
            (glass.C1 * λ^2 / (λ^2 - glass.C2)) +
            (glass.C3 * λ^2 / (λ^2 - glass.C4)) +
            (glass.C5 * λ^2 / (λ^2 - glass.C6)) +
            (glass.C7 * λ^2 / (λ^2 - glass.C8))
        n_rel = sqrt(formula_rhs + one(T))
    elseif glass.dispform == 7
        # HandbookOfOptics1
        formula_rhs =
            T(glass.C1) + (glass.C2 / (λ^2 - glass.C3)) - (glass.C4 * λ^2)
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 8
        # HandbookOfOptics2
        formula_rhs =
            T(glass.C1) + (glass.C2 * λ^2 / (λ^2 - glass.C3)) - (glass.C4 * λ^2)
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 9
        # Sellmeier4
        formula_rhs =
            T(glass.C1) +
            (glass.C2 * λ^2 / (λ^2 - glass.C3)) +
            (glass.C4 * λ^2 / (λ^2 - glass.C5))
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 10
        # Extended1
        formula_rhs =
            T(glass.C1) +
            (glass.C2 * λ^2) +
            (glass.C3 * λ^(-2)) +
            (glass.C4 * λ^(-4)) +
            (glass.C5 * λ^(-6)) +
            (glass.C6 * λ^(-8)) +
            (glass.C7 * λ^(-10)) +
            (glass.C8 * λ^(-12))
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 11
        # Sellmeier5
        formula_rhs =
            (glass.C1 * λ^2 / (λ^2 - glass.C2)) +
            (glass.C3 * λ^2 / (λ^2 - glass.C4)) +
            (glass.C5 * λ^2 / (λ^2 - glass.C6)) +
            (glass.C7 * λ^2 / (λ^2 - glass.C8)) +
            (glass.C9 * λ^2 / (λ^2 - glass.C10))
        n_rel = sqrt(formula_rhs + one(T))
    elseif glass.dispform == 12
        # Extended2
        formula_rhs =
            T(glass.C1) +
            (glass.C2 * λ^2) +
            (glass.C3 * λ^(-2)) +
            (glass.C4 * λ^(-4)) +
            (glass.C5 * λ^(-6)) +
            (glass.C6 * λ^(-8)) +
            (glass.C7 * λ^4) +
            (glass.C8 * λ^6)
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 13
        # Extended3
        formula_rhs =
            T(glass.C1) +
            (glass.C2 * λ^2) +
            (glass.C3 * λ^(4)) +
            (glass.C4 * λ^(-2)) +
            (glass.C5 * λ^(-4)) +
            (glass.C6 * λ^(-6)) +
            (glass.C7 * λ^(-8)) +
            (glass.C8 * λ^(-10)) +
            (glass.C9 * λ^(-12))
        n_rel = sqrt(formula_rhs)
    else
        @error "Invalid glass dispersion formula"
    end

    # get the absolute index of the material
    n_abs = n_rel * n_air_at_ref

    # If "TD" is included in the glass data, then include pressure and temperature dependence of the lens
    # environment. From Schott"s technical report "TIE-19: Temperature Coefficient of the Refractive Index".
    # The above "n_rel" data are assumed to be from the reference temperature T_ref. Now we add a small change
    # delta_n to it due to a change in temperature.
    ΔT = temperature - reference_temp
    if !isnan(glass.D₀) &&
       abs(ΔT) > 0.0 &&
       (
           glass.D₀ != 0 ||
           glass.D₁ != 0 ||
           glass.D₂ != 0 ||
           glass.E₀ != 0 ||
           glass.E₁ != 0
       )
        Sₜₖ = glass.λₜₖ < 0.0 ? -one(T) : one(T)
        Δn_abs =
            ((n_rel^2 - one(T)) / (2.0 * n_rel)) * (
                glass.D₀ * ΔT +
                glass.D₁ * ΔT^2 +
                glass.D₂ * ΔT^3 +
                ((glass.E₀ * ΔT + glass.E₁ * ΔT^2) / (λ^2 - Sₜₖ * glass.λₜₖ^2))
            )
        n_abs = n_abs + Δn_abs
    end

    # make the index relative to the RI of the air at the system temperature/pressure again
    n_rel = n_abs / n_air_at_sys
    return n_rel
end

# ============================================================================ #
# ====                        Material Table                              ==== #
# ============================================================================ #

MaterialTable = Dict{Symbol, AbstractMaterial}(
    :Air => Air(),
    :PerfectAbsorber => PerfectAbsorber(),
    :NBK7 => NBK7,
)

function registerMaterial(name::Symbol, mat::T) where {T <: AbstractMaterial}
    if name in MaterialTable
        error("name '$name' already registered.")
    end
    return MaterialTable[name] = mat
end

end # module
