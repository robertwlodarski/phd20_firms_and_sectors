# Content
# 1. Wrapper function 
# 2. Experiment 

# 1. Wrapper function 
function fnRunCounterfactual(params::ModelParameters, var_name::Symbol, new_value::Float64)
    
    # A. Create a new struct with the modified field
    𝕡̃       = ModelParameters(params; NamedTuple{(var_name,)}((new_value,))...)
    𝕡ᶜᶠ     = fnSetUpParameters(𝕡̃)
    eᶜᶠ     = fnSetUpEndo(𝕡ᶜᶠ)
    
    # B. Solve the model at the new parameter value
    pᶜᶠ = fnSolveSteadyState(𝕡ᶜᶠ, eᶜᶠ) 
    
    # C. Calculate the key moments
    Mᶜᶠ  = sum(eᶜᶠ.μ⃗)
    𝔼nᶜᶠ = dot(eᶜᶠ.n⃗, eᶜᶠ.μ⃗) / Mᶜᶠ
    𝔼φᶜᶠ = dot(𝕡ᶜᶠ.φ⃗, eᶜᶠ.μ⃗) / Mᶜᶠ 
    ERᶜᶠ = eᶜᶠ.M / Mᶜᶠ
    
    # D. Return results
    return (
        Variable    = string(var_name),
        Value       = new_value,
        p̂           = pᶜᶠ,
        ER          = ERᶜᶠ,
        M̂ᵢ          = Mᶜᶠ,
        𝔼n̂          = 𝔼nᶜᶠ,
        𝔼φ          = 𝔼φᶜᶠ
    )
end

# 2. Experiment 
# A. Define the baseline for comparison
baseline = (
    Variable    = "Baseline",
    Value       = 0.0,
    p̂           = p̂,
    ER          = Endo.M / sum(Endo.μ⃗),        # Changed to Endo
    M̂ᵢ          = sum(Endo.μ⃗),                 # Changed to Endo
    𝔼n̂          = dot(Endo.n⃗, Endo.μ⃗) / sum(Endo.μ⃗), # Changed to Endo
    𝔼φ          = dot(UsedParameters.φ⃗, Endo.μ⃗) / sum(Endo.μ⃗) # Changed to UsedParameters/Endo
)

# B. Define the counterfactual scenarios
scenarios = [
    (:c,   UsedParameters.c * 1.2),  
    (:cₑ,  UsedParameters.cₑ * 1.2), 
    (:D̄,   UsedParameters.D̄ * 1.2),  
    (:φ̄,   exp(1.112))               
]

# C. Run the loop
cf_results  = [fnRunCounterfactual(UsedParameters, s[1], s[2]) for s in scenarios]
all_results = [baseline; cf_results]
df_results  = DataFrame(all_results)
df_results