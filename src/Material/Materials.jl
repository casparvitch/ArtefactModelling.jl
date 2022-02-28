
module Materials

abstract type AbstractMaterial end

abstract type AbstractTransparentMaterial end
# add absorption etc. methods, they'll be useful

struct Glass <: AbstractTransparentMaterial end

struct Air <: AbstractTransparentMaterial end

struct MaterialID end

end # module
