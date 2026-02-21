# Content 
# 1. Parameters (structure and constructor)
# 2. Endogenous variables (structure and constructor)

# 1. Parameters (structure)
@with_kw struct ModelParameters

    # A. Parameters 
    α::Float64          = 1/3   # DRS parameter 
    β::Float64          = 0.96  # Discount 
    c::Float64          = 20    # Fixed cost 
    cₑ::Float64         = 39.7  # Entry cost 
    D̄::Float64          = 100   # Demand parameter
    ε::Float64          = 1.0   # Elasticity 
    N̄::Float64          = 66.67 # Labour supply 

    # B. AR(1) parameters 
    ρ::Float64          = 0.9       # Persistence (log) 
    σ::Float64          = 0.4       # Standard deviation (log)
    φ̄::Float64          = exp(1.39) # Average productivity 

    # C. Productivity grid 
    Nᵩ::Int             = 50                # Productivity grid size 
    φ⃗::Vector{Float64}  = zeros(Nᵩ)         # Productivity grid 
    Γ::Matrix{Float64}  = zeros(Nᵩ,Nᵩ)      # Transition matrix 
    ν⃗::Vector{Float64}  = zeros(Nᵩ)         # Stationary distribution 

    # D. Algorithm settings  
    δᵛᶠⁱ::Float64       = 1e-5      # Tolerance parameter for VFI 
    𝒾̄ᵛᶠⁱ::Int           = 2000      # Maximum VFI iterations  
    p̲::Float64          = 1e-3      # Minimum admissible price 
    p̅::Float64          = 10        # Maximum admissible price
end 

# 2. Parameters (constructor)
function fnSetUpParameters(params::ModelParameters = ModelParameters())

    # A. Unpacking business 
    @unpack φ̄, ρ, σ, Nᵩ = params

    # B. Idiosyncratic productivity items → Rouwenhorst as ρ ≃ 1
    # xₜ = ρ xₜ₋₁ + ϵₜ, where xₜ = log φₜ - log φ̄
    ℳ𝒞                  = rouwenhorst(Nᵩ,ρ,σ)               

    # C. Save results 
    return reconstruct(params;
        Γ               = ℳ𝒞.p,
        φ⃗               = exp.(ℳ𝒞.state_values .+ log(φ̄)),
        ν⃗               = stationary_distributions(ℳ𝒞)[1]
    )
end 
UsedParameters = fnSetUpParameters()

# 2. Endogenous variables preallocation (structure)
@with_kw mutable struct EndogenousVariables

    # A. Value functions 
    V⃗::Vector{Float64}      # Value function
    V⃗ᶜ::Vector{Float64}     # Continuation value 
    Vᵉ::Float64             # Entry value 
    π⃗::Vector{Float64}      # Profit function 
    n⃗::Vector{Float64}      # Hiring policy 
    𝕀ᵉ::Vector{Bool}        # Entry policy 
    𝕀ᶜ::Vector{Bool}        # Continuation policy
    φ̲ᵢ::Int                 # Exit threshold index  
    μ̃⃗::Vector{Float64}      # Per entrant distribution
    Γ̃::Matrix{Float64}      # Modified transition matrix
    Q̃ˢ::Float64             # Per entrant output 
    M::Float64              # Entry mass  
    D::Float64              # Demand  
    μ⃗::Vector{Float64}      # Equilibrium distribution of firms
    y⃗::Vector{Float64}      # Output per type of firm 
end

# 2. Endogenous variables preallocation (constructor)
function fnSetUpEndo(params::ModelParameters)

    # A. Unpacking business 
    @unpack Nᵩ = params 

    # B. Preallocate values 
    V⃗       = zeros(Nᵩ) 
    V⃗ᶜ      = zeros(Nᵩ)     
    Vᵉ      = 0.0
    π⃗       = zeros(Nᵩ)
    n⃗       = zeros(Nᵩ)
    φ̲ᵢ      = 1
    𝕀ᵉ      = fill(true,Nᵩ)
    𝕀ᶜ      = fill(true,Nᵩ)
    μ̃⃗       = zeros(Nᵩ)
    Γ̃       = zeros(Nᵩ,Nᵩ)     
    Q̃ˢ      = 0.0
    M       = 0.0
    D       = 0.0
    μ⃗       = zeros(Nᵩ)
    y⃗       = zeros(Nᵩ)

    # C. Return 
    return EndogenousVariables(
        V⃗   = V⃗,
        V⃗ᶜ  = V⃗ᶜ,
        Vᵉ  = Vᵉ,
        π⃗   = π⃗, 
        n⃗   = n⃗,
        φ̲ᵢ  = φ̲ᵢ,
        𝕀ᵉ  = 𝕀ᵉ,
        𝕀ᶜ  = 𝕀ᶜ,
        μ̃⃗   = μ̃⃗,
        Γ̃   = Γ̃,
        Q̃ˢ  = Q̃ˢ,
        M   = M,
        D   = D,
        μ⃗   = μ⃗,
        y⃗   = y⃗
    )   
end 
Endo    = fnSetUpEndo(UsedParameters)