module HeuristicNLSolve

using MPI, NLsolve, Format, DelimitedFiles
import Random: shuffle!, seed!
import LinearAlgebra: Transpose

# Message integer flags
const SENDING       = 2
const RESET         = 1
const CONTINUE      = 0
const STOP          = -1

const SEARCH_PROC   = 111
const SOLVE_PROC    = 222

include("SwarmOptions.jl")
include("SolverOptions.jl")
include("Problem.jl")
include("Particle.jl")
include("Swarm.jl")
include("DMSPSO.jl")
include("SolverServer.jl")

export SwarmOptions, SolverOptions, Problem, DMSPSO, SolverServer
export optimize!, start!

end
