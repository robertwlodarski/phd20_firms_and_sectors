# Contents 
# 1. Firm optimisation 
# 2. Stationary distribution 

function fnFirmOptim!(params,endo,p,w)

    # A. Unpacking business 
    @unpack β, α, γ, φ⃗, τ⃗, f = params 

    # B. Capital and labour
    r           = 1 / β - 1 
    ratio       = γ / α * w / r
    endo.𝐤      .= (ratio^α .* r ./ ((1 .- τ⃗') .* p .* φ⃗ .* γ)).^(1 / (γ + α - 1))
    endo.𝐧      .= (ratio^(-γ) .* w ./ ((1 .- τ⃗') .* p .* φ⃗ .* α)).^(1 / (γ + α - 1))

    # C. Profits and policies 
    endo.Π      .= (1 .- τ⃗') .* p .* φ⃗ .* (endo.𝐤).^(γ) .* endo.𝐧.^(α) .- r .* endo.𝐤 .- w .* endo.𝐧 .- f
    endo.𝐞      .= (endo.Π .>= 0)
end 

# 2. Stationary distribution 
function fnStationaryDist!(params, endo, m)
    
    # A. Unpacking business 
    @unpack λ, g, φ⃗, γ, α = params 

    # B. Compute the distribution 
    endo.μ      .= m / (1 - λ) .* endo.𝐞 .* g 

    # C. Aggregation 
    endo.Kᴰ     = sum(endo.𝐤 .* endo.μ)
    endo.Nᴰ     = sum(endo.𝐧 .* endo.μ)
    endo.Y      = sum( φ⃗ .* endo.𝐤.^γ .* endo.𝐧.^α .* endo.μ)
end 