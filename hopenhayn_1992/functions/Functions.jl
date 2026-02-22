# Content 
# 1. Static policy function
# 2. Initial guess of value function 
# 3. Value function iteration 
# 4. Exit threshold 
# 5. Entry error 
# 6. Updating distributions
# 7. Aggregation business
# 8. Solve for steady state 

# 1. Static policy function  
function fnStaticPolicies!(params, endo,p)

    # A. Unpacking business
    @unpack α, c, φ⃗ = params 

    # B. Update the profit function 
    endo.π⃗          .= (1 - α) * α^(α / (1 - α)) * (p .* φ⃗).^(1 / (1 - α)) .- c 

    # C. Update the labour demand 
    endo.n⃗          .= (α * p .* φ⃗).^(1 /(1 - α))
end 

# 2. Initial guess of value function 
function fnV⁰!(params, endo)

    # A. Unpacking business
    @unpack β, Γ = params 

    # # B. Prepare the initial guess (use the expected value)
    endo.V⃗          .=  endo.π⃗ .+ β .* (Γ * endo.π⃗) 
end 

# 3. Value function iteration 
function fnVFI!(params, endo)

    # A. Unpacking business
    @unpack β, Γ, δᵛᶠⁱ, 𝒾̄ᵛᶠⁱ = params

    # B. Initial guess 
    fnV⁰!(params, endo)
    V⃗ᵒˡᵈ        = endo.V⃗
    V⃗ᶜ          = similar(V⃗ᵒˡᵈ)
    V⃗ⁿᵉʷ        = similar(V⃗ᵒˡᵈ)
    ϵᵛᶠⁱ        = 1.0
    𝒾           = 1

    # C. Start the loop 
    while (ϵᵛᶠⁱ > δᵛᶠⁱ && 𝒾 < 𝒾̄ᵛᶠⁱ)
        V⃗ᶜ          .= Γ * V⃗ᵒˡᵈ
        V⃗ⁿᵉʷ        .= endo.π⃗ .+ β .* max.(0,V⃗ᶜ)
        ϵᵛᶠⁱ        = maximum(abs.(V⃗ᵒˡᵈ .- V⃗ⁿᵉʷ))
        𝒾           += 1 
        V⃗ᵒˡᵈ        .= V⃗ⁿᵉʷ

        # Warning message  
        if 𝒾 == 𝒾̄ᵛᶠⁱ
            println("VFI did not converge, lol. Final distance: ", ϵᵛᶠⁱ)
        end 
    end 

    # D. Save the result 
    endo.V⃗          .= V⃗ⁿᵉʷ
    endo.V⃗ᶜ         .= V⃗ᶜ
end 

# 4. Exit threshold 
function fnExitThreshold!(endo,params)

    # A. Unpacking business 
    @unpack Nᵩ, φ⃗, Γ = params

    # B. Find the index the first firm that decides to stay
    endo.𝕀ᶜ.= (endo.V⃗ᶜ .>= 0.0)
    idx     = findfirst(endo.𝕀ᶜ)
    endo.φ̲ᵢ = isnothing(idx) ? Nᵩ : idx    

    # C. Modified transition operator accounting for exit
    Γ̃           = copy(Γ')
    # Note the transpose: Matrix Γ̃ = Γᵀ with exit rows zeroed.
    Γ̃[:,.!endo.𝕀ᶜ]  .= 0.0
    endo.Γ̃          = Γ̃
end 

# 5. Entry error 
function fnEntryError!(endo, params, p)

    # A. Unpacking
    @unpack ν⃗,cₑ = params 

    # B. Static setting
    fnStaticPolicies!(params, endo,p)
    fnV⁰!(params, endo)

    # C. VFI 
    fnVFI!(params, endo)
    fnExitThreshold!(endo,params)

    # D. Find Vᵉ(p) 
    endo.Vᵉ = dot(endo.V⃗, ν⃗)
    return  endo.Vᵉ - cₑ
end 

# 6. Compute the per-entrant distribution 
function fnDistributions!(endo,params)

    # A. Unpacking business 
    @unpack ν⃗ = params 
    
    # B. Compute the distribution 
    endo.μ̃⃗      = (I - endo.Γ̃) \ ν⃗
end 

# 7. Aggregation business 
function fnAggregation!(endo,params,p)

    # A. Unpacking business 
    @unpack α, φ⃗, D̄,ε    = params 

    # B. Final product demand and per-entrant output
    endo.D      = D̄ / (p^ε)
    endo.y⃗      = (α * p)^(α /(1-α)) .*  φ⃗.^(1 / (1-α))
    endo.Q̃ˢ     = dot(endo.y⃗,endo.μ̃⃗)

    # C. Mass of firms
    endo.M      = endo.D / endo.Q̃ˢ

    # Equilibrium distribution 
    endo.μ⃗      = endo.M .* endo.μ̃⃗
end 

# 8. Solve for steady state 
function fnSolveSteadyState(params::ModelParameters,endo::EndogenousVariables)
    
    # A. Unpacking business 
    @unpack p̲, p̅, φ⃗, c, cₑ = params 
    
    # B. Set up the optimisation 
    𝓅           = (p̲,p̅)
    𝒻(p)        = fnEntryError!(endo, params, p)

    # C. Solve the model 
    # println("Solving for the steady state equilibrium")
    p̂           = find_zero(𝒻,𝓅)
    # println("Equilibrium price:                     $(round(p̂,digits=2))")

    # D. Recalculate everything (do I need it?)
    ϵ̂ᵉ          = fnEntryError!(endo, params, p̂)
    # println("Final error in the entry clearing:     $(round(ϵ̂ᵉ,digits=2))")

    # E. Get distributions 
    fnDistributions!(endo,params)
    fnAggregation!(endo,params,p̂)

    # F. Print messages 
    # Variables for printing
    M̂       = endo.M        # Entrants 
    M̂ⁱ      = sum(endo.μ⃗)   # Incumbents 
    𝔼n̂      = dot(endo.n⃗,endo.μ⃗) / M̂ⁱ
    𝔼ŷ      = dot(endo.y⃗,endo.μ⃗) / M̂ⁱ
    𝔼φ̂      = dot(φ⃗,endo.μ⃗) / M̂ⁱ
    𝔼π̂      = dot(endo.π⃗,endo.μ⃗) / M̂ⁱ
    𝔼φ̲̂      = φ⃗[endo.φ̲ᵢ] / maximum(φ⃗)
    𝔼N̂      = dot(endo.n⃗,endo.μ⃗) 
    𝔼Q̂      = endo.D
    # Printout 
    # println("Mass of entrants:                      $(round(M̂,digits=2))")
    # println("Average firm size:                     $(round(𝔼n̂,digits=2))")
    # println("Average firm output:                   $(round(𝔼ŷ,digits=2))")
    # println("Average firm productivity:             $(round(𝔼φ̂,digits=2))")
    # println("Average firm profit:                   $(round(𝔼π̂,digits=2))")
    # println("Entry threshold (% max productivity):  $(round(𝔼φ̲̂,digits=2))")
    # println("Total labour demand:                   $(round(𝔼N̂,digits=2))")
    # println("Total production:                      $(round(𝔼Q̂,digits=2))")

    # G. Return 
    return p̂
end 