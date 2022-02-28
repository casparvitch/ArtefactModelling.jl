
module Vis

using GLMakie: GLMakie
using ColorTypes: RGB
using FromFile: @from
@from "../Geometry/Geometry.jl" import Geometry:
    Surface,
    CSGTree,
    CSGGenerator,
    BoundingBox,
    TriangleMesh,
    makemesh,
    makiemesh

global current_main_scene = nothing
global current_layout_scene = nothing
global current_3d_scene = nothing
global current_mode = nothing
# modes:    nothing, :default  -> Original Vis beheviour    
#           :pluto             -> support pluto notebooks 
#           :docs              -> support documenter figures 

# added the following 2 functions to allow us to hack the drawing mechanisim while in a pluto notebook
set_current_main_scene(scene) = (global current_main_scene = scene)
set_current_3d_scene(lscene) = (global current_3d_scene = lscene)

get_current_mode() = begin
    global current_mode
    return current_mode
end
set_current_mode(mode) = (global current_mode = mode)

"""
    scene(resolution = (1000, 1000))

Create a new GLMakie scene with the given resolution including control buttons.
"""
function scene(resolution = (1000, 1000))
    @assert resolution[1] > 0 && resolution[2] > 0

    scene, layout = GLMakie.layoutscene(resolution = resolution)
    global current_main_scene = scene
    global current_layout_scene = layout
    lscene =
        layout[1, 1] = GLMakie.LScene(
            scene,
            scenekw = (
                camera = GLMakie.cam3d_cad!,
                axis_type = GLMakie.axis3d!,
                raw = false,
            ),
        )
    global current_3d_scene = lscene

    # in these modes we want to skip the creation of the utility buttons as these modes are not interactive
    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return scene, lscene
    end

    threedbutton = GLMakie.Button(
        scene,
        label = "3D",
        buttoncolor = RGB(0.8, 0.8, 0.8),
        height = 40,
        width = 80,
    )
    twodxbutton = GLMakie.Button(
        scene,
        label = "2D-x",
        buttoncolor = RGB(0.8, 0.8, 0.8),
        height = 40,
        width = 80,
    )
    twodybutton = GLMakie.Button(
        scene,
        label = "2D-y",
        buttoncolor = RGB(0.8, 0.8, 0.8),
        height = 40,
        width = 80,
    )
    savebutton = GLMakie.Button(
        scene,
        label = "Screenshot",
        buttoncolor = RGB(0.8, 0.8, 0.8),
        height = 40,
        width = 160,
    )

    GLMakie.on(threedbutton.clicks) do nclicks
        make3d(lscene)
        return yield()
    end

    GLMakie.on(twodybutton.clicks) do nclicks
        make2dy(lscene)
        return yield()
    end

    GLMakie.on(twodxbutton.clicks) do nclicks
        make2dx(lscene)
        return yield()
    end

    GLMakie.on(savebutton.clicks) do nclicks
        Vis.save("screenshot.png")
        return yield()
    end

    layout[2, 1] = GLMakie.grid!(
        hcat(threedbutton, twodxbutton, twodybutton, savebutton),
        tellwidth = false,
        tellheight = true,
    )
    return scene, lscene
end

function make3d(scene::GLMakie.LScene = current_3d_scene)
    s = scene.scene
    # use 3d camera
    GLMakie.cam3d_cad!(s)
    # reset scene rotation
    s.transformation.rotation[] = GLMakie.Quaternion(0.0, 0.0, 0.0, 1.0)
    # show all the axis ticks
    s[GLMakie.OldAxis].attributes.showticks[] = (true, true, true)
    # reset tick and axis label rotation and position
    rot = (
        GLMakie.Quaternion(0.0, 0.0, -0.7071067811865476, -0.7071067811865475),
        GLMakie.Quaternion(0.0, 0.0, 1.0, 0.0),
        GLMakie.Quaternion(0.0, 0.7071067811865475, 0.7071067811865476, 0.0),
    )
    s[GLMakie.OldAxis].attributes.ticks.rotation[] = rot
    s[GLMakie.OldAxis].attributes.names.rotation[] = rot
    s[GLMakie.OldAxis].attributes.ticks.align =
        ((:left, :center), (:right, :center), (:right, :center))
    s[GLMakie.OldAxis].attributes.names.align =
        ((:left, :center), (:right, :center), (:right, :center))
    # reset scene limits to automatic
    s.limits[] = GLMakie.Makie.Automatic()
    return GLMakie.update!(s)
end

function make2dy(scene::GLMakie.LScene = current_3d_scene)
    s = scene.scene
    # use 2d camera
    GLMakie.cam2d!(s)

    scene_transform = GLMakie.qrotation(GLMakie.Vec3f(0, 1, 0), 0.5pi)
    scene_transform_inv = GLMakie.qrotation(GLMakie.Vec3f(0, 1, 0), -0.5pi)    # to use with the ticks and names

    # set rotation to look onto yz plane
    s.transformation.rotation[] = scene_transform
    # hide x ticks

    # there is a bug in GLMakie 0.14.2 which cause an exception setting the X showticks to false. 
    # we work wround it by making sure the labels we want to turn off are ortogonal to the view direction 
    # s[GLMakie.OldAxis].attributes.showticks[] = (false, true, true)
    s[GLMakie.OldAxis].attributes.showticks[] = (true, true, true)

    # set tick and axis label rotation and position
    s[GLMakie.OldAxis].attributes.ticks.rotation[] =
        (0.0, scene_transform_inv, scene_transform_inv)
    s[GLMakie.OldAxis].attributes.names.rotation[] =
        s[GLMakie.OldAxis].attributes.ticks.rotation[]
    s[GLMakie.OldAxis].attributes.ticks.align =
        ((:right, :center), (:right, :center), (:center, :right))
    s[GLMakie.OldAxis].attributes.names.align =
        ((:left, :center), (:left, :center), (:center, :left))
    # update the scene limits automatically to get true reference values
    s.limits[] = GLMakie.Makie.Automatic()
    GLMakie.update_limits!(s)
    # manually set the scene limits to draw the axes correctly
    o, w = GLMakie.origin(s.data_limits[]), GLMakie.widths(s.data_limits[])
    s.limits[] = GLMakie.FRect3D((1000.0f0, o[2], o[3]), (w[2], w[2], w[3]))
    # set the eye (i.e. light) position to behind the camera
    s.camera.eyeposition[] = (0, 0, -100)
    return GLMakie.update!(s)
end

function make2dx(scene::GLMakie.LScene = current_3d_scene)
    s = scene.scene
    # use 2d camera
    GLMakie.cam2d!(s)

    scene_transform =
        GLMakie.qrotation(GLMakie.Vec3f(0, 0, 1), 0.5pi) *
        GLMakie.qrotation(GLMakie.Vec3f(1, 0, 0), 0.5pi)
    scene_transform_inv =
        GLMakie.qrotation(GLMakie.Vec3f(1, 0, 0), -0.5pi) *
        GLMakie.qrotation(GLMakie.Vec3f(0, 0, 1), -0.5pi)

    # set rotation to look onto yz plane
    s.transformation.rotation[] = scene_transform
    # hide y ticks

    # there is a bug in GLMakie 0.14.2 which cause an exception setting the X showticks to false. 
    # we work wround it by making sure the labels we want to turn off are ortogonal to the view direction 
    s[GLMakie.OldAxis].attributes.showticks[] = (true, false, true)

    # set tick and axis label rotation and position
    s[GLMakie.OldAxis].attributes.ticks.rotation[] =
        (scene_transform_inv, 0.0, scene_transform_inv)
    s[GLMakie.OldAxis].attributes.names.rotation[] =
        s[GLMakie.OldAxis].attributes.ticks.rotation[]
    s[GLMakie.OldAxis].attributes.ticks.align =
        ((:right, :center), (:right, :center), (:center, :center))
    s[GLMakie.OldAxis].attributes.names.align =
        ((:left, :center), (:right, :center), (:center, :center))
    # update the scene limits automatically to get true reference values
    s.limits[] = GLMakie.Makie.Automatic()
    GLMakie.update_limits!(s)
    # manually set the scene limits to draw the axes correctly
    o, w = GLMakie.origin(s.data_limits[]), GLMakie.widths(s.data_limits[])
    s.limits[] = GLMakie.FRect3D((o[1], -1000.0f0, o[3]), (w[1], w[1], w[3]))
    # set the eye (i.e. light) position to behind the camera
    s.camera.eyeposition[] = (0, 0, -100)
    return GLMakie.update!(s)
end

"""
    draw(ob; resolution = (1000, 1000), kwargs...)

Draw an object in a new scene.
`kwargs` depends on the object type.
"""
function draw(ob; resolution = (1000, 1000), kwargs...)
    scene, lscene = Vis.scene(resolution)
    draw!(lscene, ob; kwargs...)
    display(scene)

    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return scene
    end
end

"""
    draw!([scene = currentscene], ob; kwargs...)

Draw an object in an existing scene.
`kwargs` depends on the object type.
"""
function draw!(ob; kwargs...)
    if current_3d_scene === nothing
        scene, lscene = Vis.scene()
    else
        scene = current_main_scene
        lscene = current_3d_scene
    end

    GLMakie.cameracontrols(lscene.scene).attributes.mouse_rotationspeed[] =
        0.0001f0 #doesn't seem to do anything
    draw!(lscene, ob; kwargs...)
    display(scene)

    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return scene
    end
end

## GEOMETRY

"""
    draw!(scene::GLMakie.LScene, surf::Surface{T}; numdivisions = 20, normals = false, normalcolor = :blue, kwargs...)

Transforms `surf` into a mesh using [`makemesh`](@ref) and draws the result.
`normals` of the surface can be drawn at evenly sampled points with provided `normalcolor`.
`numdivisions` determines the resolution with which the mesh is triangulated.
`kwargs` is passed on to the [`TriangleMesh`](@ref) drawing function.
"""
function draw!(
    scene::GLMakie.LScene,
    surf::Surface{T};
    numdivisions::Int = 30,
    normals::Bool = false,
    normalcolor = :blue,
    kwargs...,
) where {T <: Real}
    mesh = makemesh(surf, numdivisions)
    if nothing === mesh
        return
    end
    draw!(scene, mesh; kwargs..., normals = false)
    if normals
        ndirs = GLMakie.Point3f.(samplesurface(surf, normal, numdivisions ÷ 10))
        norigins =
            GLMakie.Point3f.(samplesurface(surf, point, numdivisions ÷ 10))
        GLMakie.arrows!(
            scene,
            norigins,
            ndirs,
            arrowsize = 0.2,
            arrowcolor = normalcolor,
            linecolor = normalcolor,
            linewidth = 2,
        )
    end
end

"""
    draw!(scene::GLMakie.LScene, tmesh::TriangleMesh{T}; linewidth = 3, shaded = true, wireframe = false, color = :orange, normals = false, normalcolor = :blue, transparency = false, kwargs...)

Draw a [`TriangleMesh`](@ref), optionially with a visible `wireframe`. `kwargs` are passed on to [`GLMakie.mesh`](http://makie.juliaplots.org/stable/plotting_functions.html#mesh).
"""
function draw!(
    scene::GLMakie.LScene,
    tmesh::TriangleMesh{T};
    linewidth = 3,
    shaded::Bool = true,
    wireframe::Bool = false,
    color = :orange,
    normals::Bool = false,
    normalcolor = :blue,
    transparency::Bool = false,
    kwargs...,
) where {T <: Real}
    points, indices = makiemesh(tmesh)
    if length(points) > 0 && length(indices) > 0
        GLMakie.mesh!(
            scene,
            points,
            indices;
            kwargs...,
            color = color,
            shading = shaded,
            transparency = transparency,
            visible = shaded,
        )
        if wireframe
            mesh = scene.scene[end][1]
            if shaded
                GLMakie.wireframe!(
                    scene,
                    mesh,
                    color = (:black, 0.1),
                    linewidth = linewidth,
                )
            else
                GLMakie.wireframe!(
                    scene,
                    mesh,
                    color = color,
                    linewidth = linewidth,
                )
            end
        end
    end
    if normals
        @warn "Normals being drawn from triangulated mesh, precision may be low"
        norigins =
            [GLMakie.Point3f(centroid(t)) for t in tmesh.triangles[1:10:end]]
        ndirs = [GLMakie.Point3f(normal(t)) for t in tmesh.triangles[1:10:end]]
        if length(norigins) > 0
            GLMakie.arrows!(
                scene,
                norigins,
                ndirs,
                arrowsize = 0.2,
                arrowcolor = normalcolor,
                linecolor = normalcolor,
                linewidth = 2,
            )
        end
    end
end

"""
    draw!(scene::GLMakie.LScene, meshes::Vararg{S}; colors::Bool = false, kwargs...) where {T<:Real,S<:Union{TriangleMesh{T},Surface{T}}}

Draw a series of [`TriangleMesh`](@ref) or [`Surface`](@ref) objects, if `colors` is true then each mesh will be colored automatically with a diverse series of colors.
`kwargs` are is passed on to the drawing function for each element.
"""
function draw!(
    scene::GLMakie.LScene,
    meshes::Vararg{S};
    colors::Bool = false,
    kwargs...,
) where {T <: Real, S <: Union{TriangleMesh{T}, Surface{T}}}
    for i in 1:length(meshes)
        if colors
            col = indexedcolor2(i)
        else
            col = :orange
        end
        draw!(scene, meshes[i]; kwargs..., color = col)
    end
end

function draw(
    meshes::Vararg{S};
    kwargs...,
) where {T <: Real, S <: Union{TriangleMesh{T}, Surface{T}}}
    scene, lscene = Vis.scene()
    draw!(lscene, meshes...; kwargs...)
    GLMakie.display(scene)
    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return scene
    end
end

"""
    draw!(scene::GLMakie.LScene, csg::Union{CSGTree,CSGGenerator}; numdivisions::Int = 20, kwargs...)

Convert a CSG object ([`CSGTree`](@ref) or [`CSGGenerator`](@ref)) to a mesh using [`makemesh`](@ref) with resolution set by `numdivisions` and draw the resulting [`TriangleMesh`](@ref).
"""
function draw!(
    scene::GLMakie.LScene,
    csg::CSGTree{T};
    numdivisions::Int = 30,
    kwargs...,
) where {T <: Real}
    return draw!(scene, makemesh(csg, numdivisions); kwargs...)
end
function draw!(
    scene::GLMakie.LScene,
    csg::CSGGenerator{T};
    kwargs...,
) where {T <: Real}
    return draw!(scene, csg(); kwargs...)
end

"""
    draw!(scene::GLMakie.LScene, bbox::BoundingBox{T}; kwargs...)

Draw a [`BoundingBox`](@ref) as a wireframe, ie series of lines.
"""
function draw!(
    scene::GLMakie.LScene,
    bbox::BoundingBox{T};
    kwargs...,
) where {T <: Real}
    p1 = SVector{3, T}(bbox.xmin, bbox.ymin, bbox.zmin)
    p2 = SVector{3, T}(bbox.xmin, bbox.ymax, bbox.zmin)
    p3 = SVector{3, T}(bbox.xmin, bbox.ymax, bbox.zmax)
    p4 = SVector{3, T}(bbox.xmin, bbox.ymin, bbox.zmax)
    p5 = SVector{3, T}(bbox.xmax, bbox.ymin, bbox.zmin)
    p6 = SVector{3, T}(bbox.xmax, bbox.ymax, bbox.zmin)
    p7 = SVector{3, T}(bbox.xmax, bbox.ymax, bbox.zmax)
    p8 = SVector{3, T}(bbox.xmax, bbox.ymin, bbox.zmax)
    return GLMakie.linesegments!(
        scene,
        [
            p1,
            p2,
            p2,
            p3,
            p3,
            p4,
            p4,
            p1,
            p1,
            p5,
            p2,
            p6,
            p3,
            p7,
            p4,
            p8,
            p5,
            p6,
            p6,
            p7,
            p7,
            p8,
            p8,
            p5,
        ];
        kwargs...,
    )
end

"""
    save(path::String)

Save the current GLMakie scene to an image file.
"""
function save(path::String)
    # save closes the window so just display it again as a work-around
    # for some reason the size isn't maintained automatically so we just reset it manually
    size = GLMakie.size(current_main_scene)
    GLMakie.save(path, current_3d_scene.scene)
    GLMakie.resize!(current_main_scene, size)
    return display(current_main_scene)
end
function save(::Nothing) end

end # module
