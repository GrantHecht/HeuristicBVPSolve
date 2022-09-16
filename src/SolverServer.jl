# A queue of tasks used by the solver
struct TaskQueue
    taskVec::Vector{Task}

    TaskQueue() = new(Vector{Task}(undef, 0))
end

struct SolverServer{fType,nleqsType,S}
    # NL Problem struct
    prob::Problem{fType, nleqsType, S}

    # Task queue
    tasks::TaskQueue

    SolverServer(prob::Problem{fType,nleqsType,S}) where {fType,nleqsType,S} =   
        new{fType,nleqsType,S}(prob, TaskQueue())
end

# Task queue methods
Base.push!(tq::TaskQueue, t::Task)   = push!(tq.taskVec, t)
Base.pop!(tq::TaskQueue)             = pop!(tq.taskVec) 
Base.deleteat!(tq::TaskQueue, i)     = deleteat!(tq.taskVec, i)
Base.eachindex(tq::TaskQueue)        = eachindex(tq.taskVec)
Base.getindex(tq::TaskQueue, i)      = tq.taskVec[i]

# Removed finished tasks from queue
function cleanUp!(tq::TaskQueue)
    for i in reverse(eachindex(tq))
        if istaskdone(tq[i]) == true
            deleteat!(tq, i)
        end
    end
    return nothing
end

# Wait on all tasks to finish running
function wait!(tq::TaskQueue)
    for i in reverse(eachindex(tq))
        wait(tq.taskVec[i])
        pop!(tq)
    end
    return nothing
end

# Function to startup solver
function start!(s::SolverServer, opts::SolverOptions)
    # Allocate buffer
    buffer  = zeros(length(s.prob.LB))

    # Begin server loop
    done            = false
    stopMsgRecved   = false
    while !done 
        # Request message from search algorithm
        msg::Int = MPI.Recv(Int, MPI.COMM_WORLD)

        # Is message is not DONE, get new solutions
        if msg != STOP
            # Loop to recieve all messages
            for i in 1:msg
                MPI.Recv!(buffer, MPI.COMM_WORLD)

                # Spawn new task and place in queue
                push!(s.tasks, Threads.@spawn(solve!(deepcopy(buffer), s.prob.nleqs!, opts)))
            end
        else
            stopMsgRecved = true
        end

        # Cleanup task queue
        cleanUp!(s.tasks)

        # If stop message has been recieved, wait for all tasks to finish
        # and then exit loop
        if stopMsgRecved
            wait!(s.tasks)
            done = true
        end
    end
    return nothing
end

function solve!(x0, nleqs!::Function, opts::SolverOptions)
    sol = nlsolve(only_fj!(nleqs!), x0; 
        show_trace = false, factor = 3.0, ftol = 1e-8) 
    
    # Write solution to file
    if opts.solOutputFlag == true
        f = open(opts.solOutputFile, "a")
        println(f, "Success: $(sol.f_converged), Initial Guess: $x0")
        writedlm(f, Transpose(sol.zero), ",")
        close(f)
    end
    return nothing
end
