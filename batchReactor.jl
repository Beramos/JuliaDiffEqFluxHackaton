using DifferentialEquations, Plots

"""
Monod equation

``\\mu\\left(S\\right) = \\mu_{max} \\frac{K_s + S}{S}``

with ``\\mu`` the growth rate, ``S`` the substrate concentration, ``K_s`` the
half-velocity constant
"""
μ(S) = p[1] * S / (p[3] + S)

"""
solveAndPlot

This functions takes a list of parameters p, and solves a predefined ODE problem
(monod) for initial conditions u0, time tspan.
Afterwards the results are plotted.

For usage in the interactive macro
"""
function solveAndPlot(p)
    prob = ODEProblem(monod,u0,tspan,p)
    sol = solve(prob,Tsit5())
    plot(sol,
         xlabel="Concentration",
         ylabel = "Time",
         labels = ["S", "X"])
end
"""
monod

This function defines the Monod ODE.
input:
------
du: list of derivates
u: list of variables, in this case S and X
p: list of parameters, in this case μ_max, Y and Ks
t: time

the derivates are calculated using the function μ
"""
function monod(du,u,p,t)
    S,X = u
    μ_max, Y, Ks = p

    du[1] = -μ(S)*X/Y  # substrate
    du[2] = μ(S)*X  # concentration MO
end

# define the initial values of the variables S and X
# u0 = [S0, X0]
u0 = [1.0 1.0]
# define begin and end time as a tuple
# tspan = (t0, tend)
tspan = (0, 10.0)
# define the model parameters
# p = [μ_max, Y, Ks]
p = [0.5, 0.5, 1.0]

# solve and plot the ODE problem
prob = ODEProblem(monod,u0,tspan,p)
sol = solve(prob,Tsit5())
plot(sol,
     xlabel="Concentration",
     ylabel = "Time",
     labels = ["S", "X"])


"""
Simple CSTR with monod kinetic

input:
------
du: list of derivates
u: list of variables, in this case S, X and P
    S: substrate concentration
    X: biomass concentration
    P: product concentration
p: list of parameters, in this case μ_max, Y, Ks, Cp, Sin, Xin, Pin, Q
    μ_max: maximum growth rate
    Yx: fraction of substrate for biomass growth
    Yp: fraction of substrate for substrate production
    Ks: half-velocity constant
    Sin: incoming substrate concentration
    Xin: incoming biomass concentration
    Pin: incoming product concentration
    Q_out: flow in and out the CSTR
    Qx: biomass flow that is recirculated from Q_out. Inflow Q_in is equal to
    Q_in = Q_out - Qx
t: time
"""
function monodCSTR!(du,u,p,t)
    S, X, P = u
    μ_max, Yx, Yp, Ks, Sin,
        Xin, Pin, Q_out, Qx = p # define the model parameters

    du[1] = Sin*(Q_out-Qx) - S*Q_out - μ(S)*X/Yx  # substrate
    du[2] = X*Qx - X*Q_out + μ(S)*X               # concentration MO
    du[3] = Pin*(Q_out-Qx) - P*Q_out + μ(S)*X/Yp  # product
end

# define the initial values of the variables S, X, P
# u0 = [S0, X0, P0]
u0 = [1.0 1.0 0.0]
# define begin and end time as a tuple
# tspan = (t0, tend)
tspan = (0, 20.0)
# define the model parameters
p = [
    2.0, # μ_max
    0.4, # Yx
    0.6, # Yp
    1.0, # Ks
    0.5, # Sin
    1.0, # Xin
    0.0, # Pin
    1.0, # Q_out
    0.0 # Qx
]

# solve and plot the ODE problem
prob = ODEProblem(monodCSTR!,u0,tspan,p)
sol = solve(prob,Tsit5())
plot(sol,
     xlabel="Time",
     ylabel = "Concentration",
     labels = ["S", "X", "P"])
"""
wrapping the ODE solver
solve the ODE problem for a time Δt, with initial condition u0 and parameters p
the output for the last timestep is returned (row vector)
"""
function timeWrapperModel(Δt,u0,p)
    tspan = (0, Δt)
    prob = ODEProblem(monodCSTR!,u0,tspan,p)
    sol = solve(prob,Tsit5(),dense=true)
    return sol.u[end]
end

function controller_prop(e,k)
    # clipping for values below 1e-3 to avoid unrealistic values of the control
    # action, e.g. negative values
    return max(e.*k,1e-3)
end

function solve_control_loop(t_array)
    for i in 1:length(t_array)
        global u
        global output
        u = timeWrapperModel(Δt,u,p)
        p[5] = controller_prop(Psp-u[3],4.0)
        output = [output; u]
    end
    return output
end

### ------------- Step change in Qout ------------------ ###
# define the model parameters
p = [
    2.0, # μ_max
    0.4, # Yx
    0.6, # Yp
    1.0, # Ks
    0.5, # Sin
    1.0, # Xin
    0.0, # Pin
    1.0, # Q_out
    0.0 # Qx
]
# define time step for which a control action can be applied
Δt = 1e-2
tend = 10.0
t_array = 1:Δt:tend
u = [1.0 1.0 0.0]
output = u
# for every timestep ti in t_array:
# if ti is between 5 and 10, return 2, else 1
QArray  = [5.0 ≤ ti ≤ 10.0 ? 2.0 : 1.0 for ti in t_array];

# solve the ODE for t_array in increments of Δt
for i in 1:length(t_array)
    global u
    global output
    p[end-1] = QArray[i]
    u = timeWrapperModel(Δt,u,p)
    output = [output; u]
end

plot(t_array,output[2:end,:],
     lw = 2.0,
     xlabel="Time",
     ylabel = "Concentration",
     labels = ["S", "X", "P"])

### ------------- Proportional controller ------------------ ###
# define the model parameters
p = [
    2.0, # μ_max
    0.4, # Yx
    0.6, # Yp
    1.0, # Ks
    0.5, # Sin
    1.0, # Xin
    0.0, # Pin
    1.0, # Q_out
    0.0 # Qx
]
Δt = 1e-3
tend = 10.0
t_array = 1:Δt:tend
u = [1.0 1.0 0.0]
output = u

# for every timestep ti in t_array:
# if ti is between 10 and 15, return 2, else 1
Psp  = 2.0; # setpoint for S

output = solve_control_loop(t_array)

plot(t_array,output[2:end,:],
     lw = 2.0,
     xlabel="Time",
     ylabel = "Concentration",
     labels = ["S", "X", "P"]
)


# --- Calibreren met DiffEqFlux --- #
# data is één volledige simulatie


using Flux, DiffEqFlux

"""
wrapper for monodCSTR!
fixing a couple of parameters that not need to be calibrated

p: param-object
"""
function monodCSTRWrapper!(du,u,p,t)
    S, X, P = u
    p_fixed = [
            2.0, # μ_max
            0.4, # Yx
            0.6, # Yp
            1.0, # Ks
            1.0, # Xin
            0.0, # Pin
            1.0 # Q_out
        ]
    μ_max, Yx, Yp, Ks,
        Xin, Pin, Q_out, = p_fixed
    μ(S) = p_fixed[1] * S / (p_fixed[4] + S)
    Sin, Qx = p

    du[1] = Sin*(Q_out-Qx) - S*Q_out - μ(S)*X/Yx  # substrate
    du[2] = X*Qx - X*Q_out + μ(S)*X               # concentration MO
    du[3] = Pin*(Q_out-Qx) - P*Q_out + μ(S)*X/Yp  # product
end

u0 = [1.0 1.0 0.0]
p = [0.5, 1.0]
prob = ODEProblem(monodCSTRWrapper!,u0,(0,10.0),p)


# initial parameter vector
p = param([0.5, 1.0])
#params = Flux.Params([p])

function predict_rd() # our 1-layer neural network
    Array(diffeq_rd(p,prob,Tsit5(),saveat=0.1))
end

"""define the loss function
predict_rd() returns a simulation of the ODE
x[1,3] is the product concentration
2 is the setpoint, we want to have a product conentration of 2
sum(abs2, ...)  is the sum of the squared absolute difference
"""
loss_rd() = sum((predict_rd()[1,3,:].-2.0).^2)

data = Iterators.repeated((), 50)
opt = ADAM(0.1)
# log the input parameters
p_log = []
cb = function () #callback function to observe training
  display(loss_rd())
  # using `remake` to re-create our `prob` with current parameters `p`
  display(plot(solve(remake(prob,p=Flux.data(p)),Tsit5(),saveat=0.1),ylim=(0,6)))
  println(Flux.data(p))
# global p_log
  # push!(p_log, Flux.data(p))
  # display(plot!(p_log, axis=:right))
end
# Display the ODE with the initial parameter values.
cb()
Flux.train!(loss_rd, [p], data, opt, cb = cb)



"""
get the data for the problem fitting the time series of the product
concentration
"""
function get_data()
    # define the initial values of the variables S, X, P
    # u0 = [S0, X0, P0]
    u0 = [1.0 1.0 0.0]
    # define begin and end time as a tuple
    # tspan = (t0, tend)
    tspan = (0, 20.0)
    # define the model parameters
    p = [
        2.0, # μ_max
        0.4, # Yx
        0.6, # Yp
        1.0, # Ks
        0.5, # Sin
        1.0, # Xin
        0.0, # Pin
        1.0, # Q_out
        0.0 # Qx
    ]

    # solve and plot the ODE problem
    prob = ODEProblem(monodCSTR!,u0,tspan,p)
    # saving data at every 0.1 seconds
    sol = solve(prob,Tsit5(),saveat=0.1)
    return [sol.u[i][2] for i in 1:length(sol.u)]
end

data_sol = get_data()
