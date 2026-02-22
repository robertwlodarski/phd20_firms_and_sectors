# Code to generate printout of the results 

function fnPrintResults(params::ModelParameters,endo::EndogenousVariables,p̂)

    # A. Unpacking business
    @unpack φ⃗, c, cₑ,Nᵩ = params 

    # B. Prepare the results
    Mᵢ              = sum(endo.μ⃗)
    𝔼φ              = dot(φ⃗,endo.μ⃗) / Mᵢ
    𝔼n              = dot(endo.n⃗.+ c,endo.μ⃗) / Mᵢ
    Y               = dot(endo.y⃗,endo.μ⃗)
    N               = dot(endo.n⃗ .+ c,endo.μ⃗) + endo.M * cₑ
    # Top 10% earnings share 
    mᵗ              = 0.10 * Mᵢ
    mᵃ              = 0.0
    nᵃ              = 0.0
    n⃗ᵉ              = (endo.n⃗ .+ c) .* endo.μ⃗
    for i in Nᵩ:-1:1
        if mᵃ + endo.μ⃗[i] <= mᵗ
            mᵃ      += endo.μ⃗[i]
            nᵃ      += n⃗ᵉ[i]
        else 
            𝓇       = mᵗ-mᵃ
            nᵃ      += (𝓇 / endo.μ⃗[i]) * n⃗ᵉ[i]
            break
        end 
    end 

    # C. Prepare the data frame
    ResultsDF       = DataFrame(
        Variable    = [
                    "Output price (p)",
                    "Entry mass (Mₑ)",
                    "Exit threshold (φ̲)",
                    "Total mass of incumbents (Mᵢ)",
                    "Average productivity across firms (𝔼φ)",
                    "Average employment per firm (𝔼n)",
                    "Total aggregate output (Y)",
                    "Labour productivity (Y/N)",
                    "Exit rate (Mₑ/Mᵢ)",
                    "Employment share of top 10%"
                ],
        Value       = [
                    round(p̂, digits = 4),
                    round(endo.M,digits = 4),
                    round(φ⃗[endo.φ̲ᵢ], digits = 4),
                    round(Mᵢ, digits = 4),
                    round(𝔼φ, digits = 4),
                    round(𝔼n, digits = 4),
                    round(Y, digits = 4),
                    round(Y/N, digits = 4),
                    round(endo.M/Mᵢ, digits = 4),
                    round(100 * nᵃ / (sum(n⃗ᵉ)), digits = 4)
        ]
    )

    # D. Print the results 
    display(ResultsDF)
end 