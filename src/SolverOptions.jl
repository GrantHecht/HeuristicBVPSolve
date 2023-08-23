
struct SolverOptions
    fTol::Float64   # Function convergence criteria
    factor::Float64 # Factor used by NLsolve to adjust trust region

    solOutputFlag::Bool
    solOutputFile::String 

    maxWaitTime::Float64  # Max wait time in hours (wait time = time after MSPSO has stopped)

    printDebugInfo::Bool
end

function SolverOptions(;fTol = 1e-8, factor = 1.0, solOutFile = "NoOutput", maxWaitTime = 1.0, printDebugInfo = false)
    if solOutFile == "NoOutput"
        solOutputFlag = false
    else
        solOutputFlag = true
        close(open(solOutFile, "w"))
    end
    SolverOptions(fTol,factor,solOutputFlag,solOutFile,maxWaitTime,printDebugInfo)
end
