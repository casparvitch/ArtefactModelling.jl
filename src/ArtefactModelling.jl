# keep it 2D for now aye

module ArtefactModelling

# ============================================================================ #
# ====                            Imports                                 ==== #
# ============================================================================ #

using GLMakie: GLMakie
using FromFile: @from
@from "Geometry/Geometry.jl" import Geometry
@from "Vis/Vis.jl" import Vis

#initialize these caches here so they will get the correct number of threads
# from the load time environment, rather than the precompile environment. 
# The latter happens if the initialization happens in the const definition. 
# If the precompile and load environments have different numbers of threads 
# this will cause an error.
function __init__()

    # this call is to try and keep the original behevior of Makie's default 
    # backend after adding the WGLMakie backend to the package
    try
        GLMakie.activate!()
    catch e
        @warn "Unable to activate! the GLMakie backend\n$e"
    end

    for _ in 1:Threads.nthreads()
        push!(
            Geometry.threadedtrianglepool,
            Dict{DataType, Geometry.TrianglePool}((
                Float64 => Geometry.TrianglePool{Float64}()
            )),
        )
        push!(
            Geometry.threadedintervalpool,
            Dict{DataType, Geometry.IntervalPool}((
                Float64 => Geometry.IntervalPool{Float64}()
            )),
        )
    end
end

# ============================================================================ #
# ====                            Procedures                              ==== #
# ============================================================================ #

function geo_test()
    ray = Geometry.Ray([0.0, 0.0, 0.0], [0.0, 1.0, 2.0])
    for obj in [
        #         Geometry.Triangle([2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]),
        Geometry.Sphere(5.0),
        Geometry.Cylinder(5.0),
        #         Geometry.SphericalCap(5.0, π / 2.0),
        #         Geometry.Rectangle(5.0, 5.0),
    ]
        println(Geometry.surfaceintersection(obj, ray))
    end
    Vis.draw(Geometry.Sphere(5.0))
    return nothing
end

function csg_test()
    # test using ∩ etc.
    gen = Geometry.Sphere(3.0) ∪ Geometry.Cylinder(2.0, 10.0)
    #         gen = Geometry.Sphere(5.0) ∩ Geometry.Cylinder(2.0, 10.0)
    #     gen = Geometry.Sphere(5.0) - Geometry.Cylinder(2.0, 10.0)
    obj = gen() # can add transform here
    ray = Geometry.Ray([0.0, 0.0, 0.0], [0.0, 1.0, 2.0])
    println(Geometry.surfaceintersection(obj, ray))
    Vis.draw(obj)
    return nothing
end

function vis_test()
    obj = Geometry.Sphere(5.0)
    return Vis.draw(obj, resolution = (300, 300))
end

function prism_test() 
    ray = Geometry.Ray([0.0, 0.0, 0.0], [0.0, 1.0, 2.0])
    for gen in [
        Geometry.Cuboid(5.0, 5.0, 5.0),
        Geometry.BoundedCylinder(2.0, 10.0),
        Geometry.HexagonalPrism(2.0, 3.0),
        Geometry.RectangularPrism(2.0, 2.0),
        Geometry.TriangularPrism(1.0),
#         Geometry.Spider(5, 2.0, 3.0)
    ]
        obj = gen()
        println(Geometry.surfaceintersection(obj, ray))
    end
    Vis.draw(Geometry.Cuboid(5.0, 5.0, 5.0))
    return nothing
end

end # module
