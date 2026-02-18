# Content 
# 1. Static policy function

# 1. Static policy function  
function fnStaticPolicies!(params, endo)

    # A. Unpacking business
    @unpack α, c, φ⃗ = params 
    @unpack p       = endo 

    # B. Update the profit function 
    endo.π⃗          .= (1 - α) * α^(α / (1 - α)) * (p .* φ⃗).^(1 / (1 - α)) .- c 

    # C. Update  the threshold
    denominator     = (1 - α)^(1 - α) * α^α * p 
    endo.φ̲          = c^(1 - α) / denominator

    # D. Update the labour demand 
    endo.n⃗          .= (α * p .* φ⃗).^(1 /(1 - α))
end 