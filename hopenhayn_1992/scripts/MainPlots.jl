# Code to generate the key plots 
# Content:
# 1. Value function vs. productivity.
# 2. Distribution plots 
# 3. Entry and stationary distributions 

# 1. Value function vs. productivity.
function fnPlotValueFunction(params::ModelParameters, endo::EndogenousVariables)
    # A. Unpacking
    @unpack φ⃗ = params 
    φ̲ = φ⃗[endo.φ̲ᵢ]

    # B. Create the base line plot
    plt = plot(φ⃗, endo.V⃗, 
        linewidth = 2.5, 
        color = :blue,
        label = "V(φ,p) ",
        xlabel = "Productivity (φ) [Log scale]", 
        ylabel = "Firm value",
        title = "Firm value and the endogenous exit decision",
        xscale = :log10,
        legend = :topleft,
        grid = true
    )

    # C. Add the shaded "Exit Region" requested in the assignment
    vspan!(plt, [φ⃗[1], φ̲], 
        color = :red, 
        alpha = 0.15, 
        label = "Vᶜ(φ,p) < 0"
    )

    # D. Add the vertical dashed line for the exact threshold
    vline!(plt, [φ̲], 
        linestyle = :dash, 
        linewidth = 2,
        color = :red, 
        label = "φ̲ = $(round(φ̲, digits=2))"
    )
    
    # E. Add a zero-line for reference
    hline!(plt, [0.0], color=:black, linewidth=1, linestyle=:dot, label="")

    # F. Display the plot
    display(plt)
end

# 2. Distribution plots 
function fnPlotEmploymentVsFirm(params::ModelParameters, endo::EndogenousVariables)
    
    # A. Unpacking
    @unpack Nᵩ, φ⃗, c = params 

    # B. Compute bin widths for true densities on a non-uniform grid
    Δφ          = zeros(Nᵩ)
    Δφ[1]       = φ⃗[2] - φ⃗[1]
    Δφ[end]     = φ⃗[end] - φ⃗[end-1]
    for i in 2:(Nᵩ-1)
        Δφ[i]   = (φ⃗[i+1] - φ⃗[i-1]) / 2.0
    end

    # C. Compute normalized probability masses
    μ_norm      = endo.μ⃗ ./ sum(endo.μ⃗)
    n_dist      = (endo.n⃗ .+ c) .* endo.μ⃗
    emp_norm    = n_dist ./ sum(n_dist)

    # D. Convert masses to true densities
    μ_dens      = μ_norm ./ Δφ
    emp_dens    = emp_norm ./ Δφ

    # E. Create panel 1: Firm distribution
    p1 = plot(φ⃗, μ_dens, 
        linewidth   = 2.5, 
        color       = :blue, 
        label       = "Firm",
        xlabel      = "Productivity (φ) [log scale]", 
        ylabel      = "Density",
        title       = "Firm",
        xscale      = :log10,
        legend      = :topright,
        grid        = true
    )
    
    # F. Create panel 2: Employment distribution
    p2 = plot(φ⃗, emp_dens, 
        linewidth   = 2.5, 
        color       = :red, 
        label       = "Employment",
        xlabel      = "Productivity (φ) [log scale]", 
        ylabel      = "Density",
        title       = "Employment",
        xscale      = :log10,
        legend      = :topright,
        grid        = true
    )
    
    # Overlay the firm distribution on the right panel for easy comparison
    plot!(p2, φ⃗, μ_dens,
        linewidth   = 2.5,
        color       = :blue,
        linestyle   = :dash,
        label       = "Firm"
    )

    # G. Combine into a side-by-side layout
    plt = plot(p1, p2, 
        layout          = (1, 2), 
        size            = (900, 400),
        bottom_margin   = 5Plots.mm,
        left_margin     = 5Plots.mm
    )

    display(plt)
end

# 3. Entry and stationary distributions 
function fnPlotEntryVsStationary(params::ModelParameters, endo::EndogenousVariables)
    
    # A. Unpacking
    @unpack Nᵩ, φ⃗, ν⃗ = params 
    φ̲       = φ⃗[endo.φ̲ᵢ] # Extract the exit threshold

    # B. Compute bin widths for true densities on a non-uniform grid
    Δφ          = zeros(Nᵩ)
    Δφ[1]       = φ⃗[2] - φ⃗[1]
    Δφ[end]     = φ⃗[end] - φ⃗[end-1]
    for i in 2:(Nᵩ-1)
        Δφ[i]   = (φ⃗[i+1] - φ⃗[i-1]) / 2.0
    end

    # C. Compute normalized probability masses (integrating to 1)
    μ_norm      = endo.μ⃗ ./ sum(endo.μ⃗)
    ν_norm      = ν⃗ ./ sum(ν⃗) 

    # Calculate exact proportions below the exit threshold
    prop_exit_entry = sum(ν_norm[1:(endo.φ̲ᵢ-1)])
    prop_exit_stat  = sum(μ_norm[1:(endo.φ̲ᵢ-1)])

    # D. Convert masses to true densities for accurate plotting
    μ_dens      = μ_norm ./ Δφ
    ν_dens      = ν_norm ./ Δφ

    # E. Create the comparative figure
    plt = plot(φ⃗, ν_dens, 
        linewidth = 2.5, 
        color = :green, 
        label   = "Entry",
        xlabel  = "Productivity (φ) [log scale]", 
        ylabel  = "Density",
        title   = "Entry and stationary distributions",
        xscale  = :log10,
        legend  = :topright,
        grid = true
    )
    
    plot!(plt, φ⃗, μ_dens, 
        linewidth   = 2.5, 
        color       = :blue, 
        label = "Stationary"
    )

    # F. Add shading, threshold line, and annotations
    # Shading the exit region
    vspan!(plt, [φ⃗[1], φ̲], 
        color = :gray, 
        alpha = 0.15, 
        label = "Exit region"
    )

    # Exit threshold line
    vline!(plt, [φ̲], 
        linestyle   = :dash, 
        linewidth   = 2,
        color       = :red, 
        label       = "Exit threshold"
    )

    # Placement for text: left-aligned in the shaded region
    x_txt = φ⃗[1] * 1.15
    y_max = maximum([maximum(ν_dens), maximum(μ_dens)])
    
    annotate!(plt, x_txt, y_max * 0.90, text("Entry exit: $(round(prop_exit_entry * 100, digits=2))%", 9, :green, :left))
    annotate!(plt, x_txt, y_max * 0.82, text("Stat. exit: $(round(prop_exit_stat * 100, digits=2))%", 9, :blue, :left))

    display(plt)
end