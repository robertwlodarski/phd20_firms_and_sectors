# Content 
# 1. Initial value function guess 
# 2. Value function iteration 
# 3. Stationary distribution 
# 4. Entry error 
# 5. Goods market clearing 
# 6. Aggregation 
# 7. Solve for the steady state 

# 1. Initial value function and policies guesses 
function fnInitialGuess!(params, endo, p)

    # A. Unpacking business 
    @unpack α, φ⃗, c, τ, n⃗, β, ν⃗ = params 

    # B. Initial employment and profit guesses 
    # Choose the static optimisation solution
    endo.𝐍          .= (α .* p .* φ⃗).^(1 / (1 - α))
    endo.Π          .= p .* φ⃗ .* endo.𝐍.^(α) .- endo.𝐍 .- c .- τ .* max.(0,n⃗' .- endo.𝐍)

    # C. Initial value function(s) 
    endo.𝐕          .= endo.Π
    endo.𝐕ᵉ         .= endo.𝐕[:,1]

    # D. Policies 
    endo.𝕀ᶠ         .= (endo.𝐍 .< n⃗')
    endo.𝕀ᶜ         .= (endo.𝐕 .>= -τ .* n⃗')
    endo.𝕀ᵉ         .= (endo.𝐕ᵉ .>= 0.0)
end 

# 2. Value function iteration 
function fnVFI!(params, endo, p)

    # A. Unpacking business 
    @unpack δᵛᶠⁱ, 𝒾̄ᵛᶠⁱ, Γ, Nᵩ, Nₙ, α, φ⃗, c, τ, n⃗, β  = params 

    # B. Initial guess 
    fnInitialGuess!(params, endo, p)
    𝐕ᵒˡᵈ        = endo.𝐕
    𝐕ⁿᵉʷ        = similar(𝐕ᵒˡᵈ)
    𝔼𝐕          = similar(𝐕ᵒˡᵈ)
    𝐍           = endo.𝐍
    𝐍ᵢ          = endo.𝐍ᵢ
    εᵛᶠⁱ        = 1.0 
    𝒾           = 1

    # C. Start the loop 
    while (εᵛᶠⁱ > δᵛᶠⁱ && 𝒾 < 𝒾̄ᵛᶠⁱ)

        # C1. Update the expectation 
        mul!(𝔼𝐕,Γ, 𝐕ᵒˡᵈ)

        # C2. Open the loop 
        for i = 1:1:Nᵩ
            for j = 1:1:Nₙ
                obj         = p * φ⃗[i] .* n⃗.^(α) .- n⃗ .- c .- τ .* max.(0.0,n⃗[j] .- n⃗) .+ β .* 𝔼𝐕[i,:]
                val, idx    = findmax(obj)
                𝐕ⁿᵉʷ[i,j]   = max(val, -τ * n⃗[j])
                𝐍[i,j]      = n⃗[idx]
                𝐍ᵢ[i,j]     = idx
            end 
        end 
        # C3. Update errors 
        εᵛᶠⁱ        = maximum(abs.(𝐕ⁿᵉʷ .- 𝐕ᵒˡᵈ))
        𝒾           += 1
        𝐕ᵒˡᵈ        .= 𝐕ⁿᵉʷ

        # Warning message  
        if 𝒾 == 𝒾̄ᵛᶠⁱ
            println("VFI did not converge, lol. Final distance: ", εᵛᶠⁱ)
        end 
    end 

    # D. Save results 
    endo.𝐕      .= 𝐕ᵒˡᵈ
    endo.𝐍      .= 𝐍
    endo.𝐍ᵢ     .= 𝐍ᵢ
    endo.𝐕ᵉ     .= endo.𝐕[:,1]
    endo.𝕀ᶠ     .= (endo.𝐍 .< n⃗')
    endo.𝕀ᶜ     .= (endo.𝐕 .> -τ .* n⃗')
    endo.𝕀ᵉ     .= (endo.𝐕ᵉ .>= 0.0)
end 

# 3. Stationary distribution 
function fnStationaryDist!(params, endo)

    # A. Unpacking business 
    @unpack Nᵩ, Nₙ, ν⃗, δᵈⁱˢᵗ,Γ,𝒾̄ᵈⁱˢᵗ = params

    # B. Settings 
    ℳ           = 1.0
    μᵒˡᵈ        = zeros(Nᵩ, Nₙ)
    μⁿᵉʷ        = similar(μᵒˡᵈ)
    εᵈⁱˢᵗ       = 1.0
    𝒾           = 1

    # C. Loop 
    while (εᵈⁱˢᵗ > δᵈⁱˢᵗ && 𝒾 < 𝒾̄ᵈⁱˢᵗ)

        # C1. Inject fresh blood 
        fill!(μⁿᵉʷ, 0.0)
        μⁿᵉʷ[:,1]   .= ℳ .* ν⃗

        # C2. Open the for loops 
        for i = 1:1:Nᵩ
            for j = 1:1:Nₙ
                idx         = endo.𝐍ᵢ[i,j]
                μⁿᵉʷ[:,idx] .+= μᵒˡᵈ[i,j] .* endo.𝕀ᶜ[i,j] .* Γ[i,:]
            end 
        end 

        # C3. Check the error 
        εᵈⁱˢᵗ   = maximum(abs.(μⁿᵉʷ .- μᵒˡᵈ))
        μᵒˡᵈ    .= μⁿᵉʷ
        𝒾       += 1

        # D. Warning message 
        if 𝒾 == 5000
        println("WARNING: Distribution loop hit 5000. Mass is exploding!")
        end
    end 

    # D. Save the distribution 
    endo.μ      .= μᵒˡᵈ
end 

# 4. Entry error 
function fnEntryError!(params, endo, p)

    # A. Unpacking
    @unpack ν⃗,cₑ,β = params 
    
    # B. VFI 
    fnVFI!(params, endo, p)

    # C. Find the expected entry value 
    𝔼𝐕ᵉ   = β * dot(ν⃗,endo.𝐕ᵉ)
    return 𝔼𝐕ᵉ - cₑ
end 

# 5. Goods market clearing 
function fnGoodsClear!(params, endo,p)

    # A. Unpacking business 
    @unpack φ⃗, θ, α = params

    # B. Unscaled goods supply 
    Yᵘ          = sum(endo.μ .* (φ⃗ .* (endo.𝐍).^(α)))

    # C. Goods demand 
    endo.C      = θ / p 
    endo.Y      = endo.C

    # D. Mass of entrants and distribution of firms 
    endo.M      = endo.C / Yᵘ
    endo.μ      .= endo.μ .* endo.M
end

# 6. Aggregation 
function fnAggregation!(params, endo)

    # A. Unpacking business 
    @unpack cₑ,τ,n⃗,α = params 

    # B. Aggregate employment
    endo.Nᵈ     = sum(endo.𝐍 .* endo.μ) + endo.M * cₑ

    # C. Aggregate firing cost 
    endo.T      = sum(τ .* max.(0.0,n⃗' .- endo.𝐍) .* endo.μ)

    # D. TFP 
    endo.A      = endo.Y / (endo.Nᵈ^α)
end 

# 7. Solve for the steady state 
function fnSolveSteadyState!(params::ModelParameters, endo::EndogenousVariables)

     # A. Unpacking business 
    @unpack p̲, p̅, φ⃗, c, cₑ,n⃗ = params 
    
    # B. Set up the optimisation 
    𝓅           = (p̲,p̅)
    𝒻(p)        = fnEntryError!(params, endo, p)

    # C. Solve the model 
    p̂           = find_zero(𝒻,𝓅)

    # D. Find stationary distribution & impose goods clearing 
    fnStationaryDist!(params, endo)
    fnGoodsClear!(params, endo,p̂)
    fnAggregation!(params, endo)

    # E. Print messages 
    # Variables for printing
    M̂       = endo.M                  # Entrants 
    M̂ⁱ      = sum(endo.μ)             # Incumbents 
    𝔼n̂      = sum(n⃗' .* endo.μ) / M̂ⁱ  # Average firm size (current)
    𝔼φ̂      = sum(φ⃗ .* endo.μ) / M̂ⁱ   # Average firm productivity
    ER      = M̂ / M̂ⁱ                  # Stationary exit rate matches entry/incumbent ratio
    
    # Printout 
    println("====================================================")
    println("         HOPENHAYN-ROGERSON (1993) EQUILIBRIUM      ")
    println("====================================================")
    println("Equilibrium price (p):                     $(round(p̂, digits=4))")
    println("Mass of entrants (M):                      $(round(M̂, digits=4))")
    println("Mass of incumbents (Mᵢ):                   $(round(M̂ⁱ, digits=4))")
    println("Stationary exit rate:                      $(round(ER * 100, digits=2))%")
    println("----------------------------------------------------")
    println("Average firm size (𝔼n̂):                    $(round(𝔼n̂, digits=2))")
    println("Average firm productivity (𝔼φ̂):            $(round(𝔼φ̂, digits=2))")
    println("----------------------------------------------------")
    println("Total production (Y):                      $(round(endo.Y, digits=2))")
    println("Total labour demand (N_d):                 $(round(endo.Nᵈ, digits=2))")
    println("Total firing costs paid (T):               $(round(endo.T, digits=4))")
    println("Aggregate TFP (A):                         $(round(endo.A, digits=4))")
    println("====================================================")
    return p̂
end 