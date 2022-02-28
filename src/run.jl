
module Run
using FromFile: @from

@from "./ArtefactModelling.jl" import ArtefactModelling

# ArtefactModelling.geo_test()
# ArtefactModelling.vis_test()
# ArtefactModelling.csg_test()
ArtefactModelling.prism_test()

end # module
