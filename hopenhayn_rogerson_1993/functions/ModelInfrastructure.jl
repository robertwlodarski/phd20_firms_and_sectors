# Content 
# 1. Parameters (structure and constructor)
# 2. Endogenous variables (structure and constructor)

# 1. Parameters (structure)
@with_kw struct ModelParameters

    # A. Parameters 
    α::Float64          = 0.64  # DRS parameter 
    β::Float64          = 0.96  # Discount 
    c::Float64          = 20    # Fixed cost 
    cₑ::Float64         = 39.7  # Entry cost 
    D̄::Float64          = 100   # Demand parameter
    ε::Float64          = 1.0   # Elasticity 
    τ::Float64          = 0.2   # Firing costs 
    N̄::Float64          = 66.67 # Labour supply 
    θ::Float64          = 1.0   # Utility of consumption multiplier 

    # B.Productivity 
    ρ::Float64          = 0.9           # Persistence (log) 
    σ::Float64          = 0.2           # Standard deviation (log)
    φ̄::Float64          = exp(1.39)     # Average productivity 
    Nᵩ::Int             = 15            # Productivity grid size 
    φ⃗::Vector{Float64}  = zeros(Nᵩ)     # Productivity grid 
    Γ::Matrix{Float64}  = zeros(Nᵩ,Nᵩ)  # Transition matrix 
    ν⃗::Vector{Float64}  = zeros(Nᵩ)     # Stationary distribution 
    Nₙ::Int             = 50            # Labour grid size 
    n⃗::Vector{Float64}  = zeros(Nₙ)     # Labour grid 
    n̅::Float64          = 10.0          # Maximum labour
    n̲::Float64          = 0.001         # Minimum labour 

    # D. Algorithm settings  
    δᵛᶠⁱ::Float64       = 1e-5      # Tolerance parameter for VFI 
    𝒾̄ᵛᶠⁱ::Int           = 2000      # Maximum VFI iterations  
    p̲::Float64          = 0.5       # Minimum admissible price 
    p̅::Float64          = 4.0       # Maximum admissible price
    δᵈⁱˢᵗ::Float64      = 1e-3      # Tolerance parameter for stationary distribution
    𝒾̄ᵈⁱˢᵗ::Int          = 5000      # Maximum distribution iterations    

end 

# 2. Parameters (constructor)
function fnSetUpParameters(params::ModelParameters = ModelParameters())

    # A. Unpacking business 
    @unpack φ̄, ρ, σ, Nᵩ, Nₙ, n̅, n̲ = params

    # B. Idiosyncratic productivity items → Rouwenhorst as ρ ≃ 1
    # xₜ = ρ xₜ₋₁ + ϵₜ, where xₜ = log φₜ - log φ̄
    ℳ𝒞                  = rouwenhorst(Nᵩ,ρ,σ)               
        
    # D. Save results 
    return reconstruct(params;
        Γ               = ℳ𝒞.p,
        φ⃗               = exp.(ℳ𝒞.state_values .+ log(φ̄)),
        ν⃗               = stationary_distributions(ℳ𝒞)[1],
        n⃗               = exp.(range(log(n̲),log(n̅), length = Nₙ))
    )
end 
UsedParameters = fnSetUpParameters()

# 2. Endogenous variables preallocation (structure)
@with_kw mutable struct EndogenousVariables

    # A. Value and policy functions 
    𝐕::Matrix{Float64}      # Value function
    𝐕ᶜ::Matrix{Float64}     # Continuation value 
    𝐕ᵉ::Vector{Float64}     # Entry value 
    Π::Matrix{Float64}      # Profit function 
    𝐍::Matrix{Float64}      # Hiring policy 
    𝐍ᵢ::Matrix{Int}         # Hiring policy (index)
    𝕀ᵉ::Vector{Bool}        # Entry policy 
    𝕀ᶜ::Matrix{Bool}        # Continuation policy
    𝕀ᶠ::Matrix{Bool}        # Firing policy

    # B. Aggregates and distributions 
    Nˢ::Float64             # Labour supply 
    M::Float64              # Mass of entrants 
    Mⁱ::Float64             # Mass of incumbents  
    Nᵈ::Float64             # Labour demand 
    Y::Float64              # Total output 
    C::Float64              # Total consumption 
    T::Float64              # Total firing cost 
    A::Float64              # TFP measurement
    μ::Matrix{Float64}      # Stationary distribution 
end

# 2. Endogenous variables preallocation (constructor)
function fnSetUpEndo(params::ModelParameters)

    # A. Unpacking business 
    @unpack Nᵩ, Nₙ = params 

    # B. Preallocate values: Values and policies 
    𝐕       = zeros(Nᵩ, Nₙ) 
    𝐕ᶜ      = zeros(Nᵩ, Nₙ)     
    𝐕ᵉ      = zeros(Nᵩ)
    Π       = zeros(Nᵩ,Nₙ)
    𝐍       = zeros(Nᵩ,Nₙ)
    𝐍ᵢ      = ones(Nᵩ,Nₙ)
    𝕀ᵉ      = fill(true,Nᵩ)
    𝕀ᶜ      = fill(true,Nᵩ,Nₙ) 
    𝕀ᶠ      = fill(true,Nᵩ,Nₙ)

    # C. Aggregates and distributions  
    Nˢ      = 0.0
    M       = 0.0
    Mⁱ      = 0.0
    Nᵈ      = 0.0
    Y       = 0.0
    C       = 0.0
    T       = 0.0
    A       = 0.0
    μ       = zeros(Nᵩ,Nₙ)

    # D. Return 
    return EndogenousVariables(
        𝐕   = 𝐕,
        𝐕ᶜ  = 𝐕ᶜ,
        𝐕ᵉ  = 𝐕ᵉ,
        Π   = Π,
        𝐍   = 𝐍,
        𝕀ᵉ  = 𝕀ᵉ,
        𝕀ᶜ  = 𝕀ᶜ,
        Nˢ  = Nˢ,
        M   = M,
        Mⁱ  = Mⁱ,
        Nᵈ  = Nᵈ,
        Y   = Y,
        T   = T,
        A   = A,
        μ   = μ,
        𝕀ᶠ  = 𝕀ᶠ,
        𝐍ᵢ  = 𝐍ᵢ,
        C   = C
    )   
end 
Endo    = fnSetUpEndo(UsedParameters)