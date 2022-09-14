# Abstract MSPSO
abstract type MSPSO end

# DMSPSO - {D}istributed {M}ulti-{S}warm {P}article {S}warm {O}ptimization
mutable struct DMSPSO{T,S,fType, C} <: MSPSO
    # Optimization problem
    prob::Problem{fType, S}

    # Process rank
    rank::Int

    # Swarm 
    swarm::Swarm{T}

    # Reset counter (Only used on master proc)
    resetCount::Int

    # PSO specific parameters
    inertiaRange::Tuple{T,T}
    minNeighborFrac::T
    selfAdjustWeight::T
    socialAdjustWeight::T
    commIterationBuffer::Int

    # MPI Swarm Communicator
    swarmComm::C
end

# DMSPSO Status Packet
# A struct for storing MSPSO status info aquired
# during communication. Status into is used 
# for status updates through IO.
mutable struct StatusPacket{T}
    globalBest::T
    swarmBests::Vector{T}
    swarmResets::Vector{Int}

    function StatusPacket(numprocs)
        T   = typeof(1.0)
        new{T}(T(Inf), Vector{T}(undef, numprocs), zeros(Int, numprocs))
    end
end

# Constructor
function DMSPSO(prob::Problem{fType,S}; numParticlesPerSwarm, inertiaRange = (0.1, 1.1),
        minNeighborFrac = 0.25, selfAdjustWeight = 1.49, socialAdjustWeight = 1.49, 
        commIterationBuffer = 1, rngSeed = nothing, comm = MPI.COMM_WORLD) where {S,fType}

    # Error checking
    length(inertiaRange) == 2 || throw(ArgumentError("inertiaRange must be of length 2."))
    minNeighborFrac > 0       || throw(ArgumentError("minNeighborFrac must be > 0."))

    # Get process rank
    rank = MPI.Comm_rank(comm)

    # Set process RNG seed
    if rngSeed === nothing
        seed!(Int(floor(time())) + rank)
    else 
        seed!(rngSeed)
    end

    # Type info
    T       = typeof(1.0)
    C       = typeof(comm)
    nIRange = (T(inertiaRange[1]), T(inertiaRange[2]))

    # Initialize swarm
    N       = length(prob.LB)
    swarm   = Swarm{T}(N, numParticlesPerSwarm)

    DMSPSO{T,S,fType, C}(prob, rank, swarm, 0, nIRange, minNeighborFrac, 
        selfAdjustWeight, socialAdjustWeight, commIterationBuffer, comm)
end

# Methods
optimize!(alg, opts) = _optimize!(alg, opts)
function _optimize!(mspso::DMSPSO, opts::SwarmOptions)
    # Check communication
    checkComm(mspso, opts)

    # Initialize swarms
    initialize!(mspso, opts)

    # Perform iteration
    iterate!(mspso, opts)
end

# Function to check all communication to try and avoid deadlocks
#   Currently, this will definitely not catch all problems but
#   hopefully it will avoid some.
function checkComm(mspso::DMSPSO, opts::SwarmOptions)
    # Need to add checks but will need to think on this a bit...
    return nothing
end

function initialize!(mspso::DMSPSO, opts::SwarmOptions)
    # Initialize swarm position and velocities
    uniformInitialization!(mspso.swarm, mspso.prob, opts)

    # Evaluate objective functions
    feval!(mspso.swarm, mspso.prob.f, opts; init = true)

    # Set swarm global best
    mspso.swarm.b = Inf
    setGlobalBest!(mspso.swarm)

    # Initialize neighborhood size
    mspso.swarm.n = max(2, floor(length(mspso.swarm) * mspso.minNeighborFrac))

    # Initialize inertia
    if mspso.inertiaRange[2] > 0
        mspso.swarm.w = (mspso.inertiaRange[2] > mspso.inertiaRange[1]) ? 
                    mspso.inertiaRange[2] : mspso.inertiaRange[1]
    else
        mspso.swarm.w = (mspso.inertiaRange[2] < mspso.inertiaRange[1]) ? 
                    mspso.inertiaRange[2] : mspso.inertiaRange[1]
    end

    # Initialize self and social adjustment
    mspso.swarm.y₁ = mspso.selfAdjustWeight
    mspso.swarm.y₂ = mspso.socialAdjustWeight
    return nothing
end

function iterate!(mspso::DMSPSO, opts::SwarmOptions)
    # Initialize time and iteration counters
    total_t0    = time()
    t0          = total_t0
    stall_t0    = t0
    comm_iters  = 0
    iters       = 0
    stallIters  = 0
    hasStalled  = false
    fStall      = Inf

    # Instantiate status packet (only used by master proc)
    if mspso.rank == 0
        statusPack = StatusPacket(MPI.Comm_size(mspso.swarmComm))
    else
        statusPack = StatusPacket(0)
    end

    # Allocate buffer for communication
    N           = length(mspso.prob.LB)
    if mspso.rank == 0
        buffer = zeros(N + 2, MPI.Comm_size(mspso.swarmComm))
    else
        buffer = zeros(N + 2, 1)
    end

    # Compute minimum neighborhood size
    minNeighborSize = max(2, floor(length(mspso.swarm)*mspso.minNeighborFrac))

    # Begin loop
    exitFlag    = 0
    resetFlag   = 0
    while exitFlag == 0
        # Inner loop to perform iterations before process communication  
        for i in 1:mspso.commIterationBuffer
            # Update iteration counter
            iters += 1

            # Prepare for evaluation
            updateVelocities!(mspso.swarm)
            step!(mspso.swarm)
            if !(all(isinf.(mspso.prob.LB)) && all(isinf.(mspso.prob.UB)))
                enforceBounds!(mspso.swarm, mspso.prob.LB, mspso.prob.UB)
            end

            # Evaluate objective function
            feval!(mspso.swarm, mspso.prob.f, opts)

            # Update global best in swarm
            flag::Bool = setGlobalBest!(mspso.swarm)

            # Update stall counter and neighborhood
            if flag
                mspso.swarm.c = max(0, mspso.swarm.c - 1)
                mspso.swarm.n = minNeighborSize
            else
                mspso.swarm.c += 1
                mspso.swarm.n = min(mspso.swarm.n + minNeighborSize,
                                    length(mspso.swarm) - 1)
            end

            # Update inertia
            if mspso.swarm.c < 2
                mspso.swarm.w *= 2.0
            elseif mspso.swarm.c > 5
                mspso.swarm.w /= 2.0
            end

            # Ensure inertia is inbounds
            if mspso.swarm.w > mspso.inertiaRange[2]
                mspso.swarm.w = mspso.inertiaRange[2]
            elseif mspso.swarm.w < mspso.inertiaRange[1]
                mspso.swarm.w = mspso.inertiaRange[1]
            end

            # Track stalling
            if fStall - mspso.swarm.b > opts.funcTol 
                fStall      = mspso.swarm.b
                stallIters  = 0
                stall_t0    = time()
            else
                stallIters += 1
            end

            # Check for swarm stopping conditions
            if stallIters >= opts.maxStallIters
                resetFlag   = true
                hasStalled  = true
                stallIters  = 0
                break
            elseif time() - t0 >= opts.maxTime
                resetFlag   = true
                break
            elseif time() - stall_t0 >= opts.maxStallTime
                resetFlag   = true
                break
            end
        end

        # Perform communication     
        if mspso.rank == 0
            statusFlag = communicateMaster!(mspso, opts, resetFlag, buffer, total_t0, statusPack) 
        else
            statusFlag = communicateWorker!(mspso, opts, resetFlag, buffer) 
        end

        # Increment communication iteration counter
        comm_iters += 1

        # Print it master proc
        if mspso.rank == 0 && opts.display == true
            printStatus(statusPack, comm_iters, (time() - total_t0)/3600.0)
        end

        # Reset or stop if commanded
        if statusFlag == RESET
            reset!(mspso, opts)
            t0          = time()
            stall_t0    = t0
            iters       = 0
            stallIters  = 0
            hasStalled  = 0
            fStall      = Inf

        elseif statusFlag == STOP
            exitFlag    = 1
        end
    end
    return nothing
end

# Here buffer is an array of floats for either sending or recieving data
# If rank == 0:
#   buffer = zeros(numdims + 2, numprocs)
# If rank != 0:
#   buffer = zeros(numdims + 2, 1)
function communicateMaster!(mspso, opts, resetFlag, buffer, total_t0, statusPack::StatusPacket)
    # Grab local buffer from buffer
    localBuffer = @view(buffer[:,1])

    # Get number of processes
    numprocs    = MPI.Comm_size(mspso.swarmComm)

    # Request info from each process and place in buffer
    for i in 1:numprocs - 1
        data, status = MPI.Recv!(localBuffer, mspso.swarmComm, MPI.Status)
        buffer[:,status.tag + 1] .= localBuffer
    end
    N                = length(mspso.prob.LB)
    buffer[1:N,1]   .= mspso.swarm.d
    buffer[N + 1, 1] = mspso.swarm.b
    buffer[N + 2, 1] = (resetFlag == true ? 1.0 : 0.0)

    # Compute resets
    resets = computeResets(mspso, opts, buffer)

    # Increment reset counter
    for i in eachindex(resets)
        if resets[i] == true
            mspso.resetCount += 1
        end
    end

    # Check for stopping conditions (Total max time option not added yet)
    stop = false
    if time() - total_t0 > opts.maxTotalTime || mspso.resetCount >= opts.maxResets 
        stop = true
    end

    # Process new solutions from resets
    #   Writes new solutions to file and send them to solver 
    #   if desired
    processSolutions(mspso, opts, resets, buffer; stop = stop)

    # Update statusPack
    fbest = Inf
    for i in 1:numprocs
        if buffer[N + 1, i] < fbest
            fbest = buffer[N + 1, i]
        end
    end
    statusPack.globalBest = fbest
    for i in eachindex(resets)
        statusPack.swarmBests[i]   = buffer[N + 1, i]  
        statusPack.swarmResets[i] += (resets[i] == true ? 1 : 0)
    end

    # Send status flags
    for i in 1:numprocs - 1
        if stop == true
            MPI.Send(STOP, mspso.swarmComm; dest = i)
        elseif resets[i] == false
            MPI.Send(CONTINUE, mspso.swarmComm; dest = i)
        else
            MPI.Send(RESET, mspso.swarmComm; dest = i)
        end
    end
    
    # Set master status flag
    if stop == true
        flag = STOP
    elseif resets[1] == false
        flag = CONTINUE
    else
        flag = RESET
    end
    return flag
end

function communicateWorker!(mspso, opts, resetFlag, buffer)
    # Fill buffer
    N = length(mspso.prob.LB) 
    buffer[1:N]    .= mspso.swarm.d
    buffer[N + 1]   = mspso.swarm.b
    buffer[N + 2]   = (resetFlag == true ? 1.0 : 0.0)

    # Send buffer to master
    MPI.Send(buffer, mspso.swarmComm; dest = 0, tag = mspso.rank)

    # Recieve reset message from master
    flag::Int = MPI.Recv(Int, mspso.swarmComm; source = 0) 
    return flag
end

# This function computes the distance between each swarm and determines if the 
# swarm should be reset
function computeResets(mspso, opts, buffer)
    # Instantiate vector of bools
    numprocs = MPI.Comm_size(mspso.swarmComm)
    resets = Vector{Bool}(undef, numprocs)

    # Set resets best on criteria not involving distance
    N = length(mspso.prob.LB)
    for i in eachindex(resets)
        if buffer[N + 2, i] > 0.9
            resets[i] = true
        else
            resets[i] = false
        end
    end

    # Compute distance between each swarm and determine if reset is necessary
    for i in 1:numprocs - 1
        for j in i + 1:numprocs
            # Check that one of the swarms is not already being reset
            if resets[i] == false && resets[j] == false
                # Compute distance between i^th and j^th swarm
                dotProd = 0
                for k in 1:N; dotProd += (buffer[k,i] - buffer[k,j])^2; end
                dist    = sqrt(dotProd)

                # Check if distance is smaller than reset distance
                if dist < opts.resetDistance
                    fi = buffer[N + 1, i]
                    fj = buffer[N + 1, j]
                    if fi < fj 
                        resets[j] = true
                    else
                        resets[i] = true
                    end
                end
            end
        end
    end

    return resets
end

function processSolutions(mspso, opts, resets, buffer; stop = false)
    if stop == false
        # Check to see if we have any resets occuring
        if any(resets) 
            # Write to file if desired
            if opts.fileOutput == true
                # Open file
                f   = open(opts.solOutFile, "a")

                # Write to file
                N   = length(mspso.prob.LB)
                for i in eachindex(resets)
                    if resets[i] == true
                        writedlm(f, Transpose(@view(buffer[1:N+1,i])), ",")
                    end
                end
                close(f)
            end

            # Send new solutions to solver process
            if opts.solverComm == true
                # Send number of new solutions
                MPI.Send(sum(resets), MPI.COMM_WORLD; dest = opts.solverRank)

                # Send solutions
                for i in eachindex(resets)
                    if resets[i] == true
                        MPI.Send(@view(buffer[1:N,i]), MPI.COMM_WORLD; dest = opts.solverRank)
                    end
                end
            end
        end
    else
        # Write all final swarm solutions to file if desired
        if opts.fileOutput == true
            f = open(opts.solOutFile, "a")
            N = length(mspso.prob.LB)
            writedlm(f, Transpose(@view(buffer[1:N+1,:])), ",")
            write(f, "DONE\n")
            close(f)
        end

        # Send all final swarm solutions to solver process
        if opts.solverComm == true
            # Send the number of swarms
            MPI.Send(MPI.Comm_size(mspso.swarmComm), MPI.COMM_WORLD; dest = opts.solverRank)

            # Send each solution
            for i in 1:MPI.Comm_size(mspso.swarmComm)
                MPI.Send(@view(buffer[1:N,i]), MPI.COMM_WORLD; dest = opts.solverRank)
            end

            # Notify solver process that we're finished over here
            MPI.Send(STOP, MPI.COMM_WORLD; dest = opts.solverRank)
        end
    end
    return nothing
end

function reset!(mspso::DMSPSO, opts::SwarmOptions)
    # Reinitialize swarm 
    uniformInitialization!(mspso.swarm, mspso.prob, opts)

    # Evaluate objective functions
    feval!(mspso.swarm, mspso.prob.f, opts; init = true)

    # Set swarm global best 
    mspso.swarm.b = Inf
    setGlobalBest!(mspso.swarm)

    # Initialize neighborhood size
    mspso.swarm.n = max(2, floor(length(mspso.swarm)*mspso.minNeighborFrac))
    return nothing
end

function printStatus(statusPack::StatusPacket, iter, time)
    totalResets = sum(statusPack.swarmResets)
    fspec1  = FormatExpr("Time Elapsed: {1:e} hours, Comm. Iterations: {2:d}")

    # Construct expresion for fbest
    expr2   = "Global Best: {1:e}, Swarm Bests: ["
    for i in eachindex(statusPack.swarmBests)
        expr2 *= "{$(i+1):e}"
        if i != lastindex(statusPack.swarmBests)
            expr2 *= ", "
        end
    end
    expr2   *= "]"
    fspec2  = FormatExpr(expr2)

    # Construct expresion for resets
    expr3   = "Total Resets: {1:d}, Swarm Resets: ["
    for i in eachindex(statusPack.swarmResets)
        expr3 *= "{$(i+1):d}"
        if i != lastindex(statusPack.swarmResets)
            expr3 *= ", "
        end
    end
    expr3   *= "]"
    fspec3  = FormatExpr(expr3)

    printfmtln(fspec1, time, iter)    
    printfmtln(fspec2, statusPack.globalBest, statusPack.swarmBests...)
    printfmtln(fspec3, totalResets, statusPack.swarmResets...)
    println("")
    return nothing
end
