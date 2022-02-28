## General Notes

  - [FromFile](https://juliapackages.com/p/fromfile)
    
      + `using FromFile; @from "file2.jl" import foo`
      + only use for internal imports etc.
      + use `using X: blaa` for registered packages etc.

  - [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl)
    
      + `using JuliaFormatter; format(".", import_to_using=true)`
      + run with: `;./formatter.sh` (at root of ArtefactModelling)
