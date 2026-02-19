# Content 
# 1. Static policy function
# 2. Initial guess of value function 
# 3. Value function iteration 
# 4. Exit threshold 
# 5. Entry error 
# 6. Finding p 
# 7. Updating distributions
# 8. Aggregation business

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
            println("Warning: VFI did not converge, lol. Final distance: ", ϵᵛᶠⁱ)
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
function fnEntryError(endo, params,p)

    # A. Static setting
    fnStaticPolicies!(params, endo,p)
    fnV⁰!(params, endo)

    # B. VFI 
    fnVFI!(params, endo)
    fnExitThreshold!(endo,params)

    # C. Find Vᵉ(p) [To be continued]

end 

# 6. Find p [To be continued]

# 7. Compute the per-entrant distribution 
function fnDistributions!(endo,params)

    # A. Unpacking business 
    @unpack ν⃗ = params 
    
    # B. Compute the distribution 
    endo.μ̃⃗      = (I - endo.Γ̃) \ ν⃗
end 

# 8. Aggregation business 
function fnAggregation(endo,params,p)

    # A. Unpacking business 
    @unpack α, φ⃗    = params 

    # B. Final product supply 
end 