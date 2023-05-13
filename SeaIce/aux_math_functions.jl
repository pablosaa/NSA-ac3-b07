# *************************************************************
# ***** MATH FUNCTIONS 
"""
function returns the truncated normal distribution of X
 julia> 𝑓ₗᵤ = 𝑁ₗᵤ(x, L=0, U=100)

Where:
* x : Vector containg data to fit distribution
* L : (optional) lower bound, default = 0
* U : (optional) upper bound, default = ∞
"""
function 𝑁ₗᵤ(x::Vector; L=0, U=nothing)
    𝑁 = let F=filter(!isnan, x)
        μ, σ = length(F)≤2 ? (NaN, 0) : (mean(F), std(F))
        
        try 
                Normal(μ, σ)
        catch e
                println("F = ", F)
                println(μ, " and ", σ)
                @error "$e"
        end
    end

    if isnothing(L) && isnothing(U)
        return 𝑁
    elseif isnothing(U)
        return truncated(𝑁, lower=L)
    elseif L<U
        return truncated(𝑁, lower=L, upper=U)
    else
        @error "L must be lower than U, but given $L ≥ $U"
    end
    
end
# ----/

"""
Function returns the mean and standard deviation of a
 truncated normal distribution of given data
julia> 
"""
function stats_𝑁ₗᵤ(x::Vector; L=0, U=nothing)
    f = 𝑁ₗᵤ(x, L=L, U=U)
    #if isnothing(L) && isnothing(U)
    #    x
    #else
    #    𝑁ₗᵤ(x, L=L, U=U)
    #end
    return mean(f), std(f), quantile(f, [.05, .25, .5, .75, .95])
            
end
# ----/


"""
faktor(x,y) computes the ratio of x/y
"""
faktor(x, y) = x/y
# ----/

"""
To compute uncertainties from function f = x*y or f=x/y 
δf(x, y, δx, δy, f) with x ± δx, y ± δy and f::Function
"""
δf(x, y, δx, δy, f) = f(x,y)*sqrt( (δx/x)^2 + (δy/y)^2 )
# ----/

"""
Convert input to SIC scaled to reference SIC given by factor fₖ,
truncating the output to 100%
"""
fix(x, fₖ) = min(fₖ*x, 100)
# ----/

"""
Function search and retrieve fix factors (μ, σ) for ASI SIC
"""
function getSICfixfactor(heute::DateTime, df::DataFrame; qq=nothing)

    vars = [:time, :ratio, :ratio_error]
    
    if !isnothing(qq) && typeof(qq)<:Vector
        foreach(qq) do V
            push!(vars, V)
        end
    end

    tmp = let heute_df = filter(:time=> X-> Date(X) == Date(heute), df)
              isempty(heute_df) && @error "No SIC fix factor on $(heute)"
              idx = abs.(heute_df[!, :time] .- heute) |> sortperm 
              heute_df[idx, vars]
    end
    
    return transform!(tmp, vars.=> ByRow(x->typeof(x)<:String ? eval(Meta.parse(x)) : x) .=> vars)

    #tmp = filter(:time=> X-> round(X, Hour(6))==round(heute, Hour(6)), df) |> F->F[!, vars]
    #return transform!(tmp, X->typeof(X)<:String ? eval(Meta.parse(X)) : X, tmp)
    #    tmp[!,N] |> X->typeof(X)<:String ? Meta.parse(X) |> eval : X
    #end
    
    #return filter(:time=> D->DateTime(D)==(heute), df) |> F->F[!, [:time, :ϱ, :δϱ]]
end
# ----/

# *******************************************************************************
"""
Transform lat, lon data to coordinates at a polar steregraphic system

julia> x, y = latlon2polarstereographic(ϕ, λ, ϕc, λ₀)

Where:
* ϕ : latitude in degree North
* λ : longitude in degree East
* x : map coordiante along X-axis [m]
* y : map coordinate along Y-axis [m]
Oprional input parameters:
* ϕ_c: latitude of true scale in degree North, standard parallel (default is -70)
* λ₀: meridian in degree East along positive Y-axis of map (default is 0)
* 𝑎 : Earth radius defined in the projection (default 6378137.0 m, WGS84)
* 𝑒 : eccentricity (default 0.08181919)

Polar stereographic projection used for satellite polar sea ice studies.
Equations adapted from "Map Projections - A Working Manual! by J.P. Snyder (1987)
"""
function latlon2xy(ϕ::T, λ::T; ϕc::T=70.0, λ₀::T=0.0, 𝑎=6378137.0, 𝑒=0.08181919) where T<:Real
    ϕ *= π/180
    ϕc*= π/180
    λ *= π/180
    λ₀*= π/180
    
    t(α) = tan(π/4-α/2)/√((1-𝑒*sin(α))/(1+𝑒*sin(α)))^𝑒  
    m(α) = cos(α)/√(1-(𝑒*sin(α))^2)
    ρ = 𝑎*m(ϕc)*t(ϕ)/t(ϕc)
    x = ρ*sin(λ-λ₀)
    y = -ρ*cos(λ-λ₀)

    return x, y
end

function latlon2xy(ϕ::Array{T}, λ::Array{T}; ϕc::T=70.0, λ₀::T=0.0, 𝑎=6378137.0, 𝑒=0.08181919) where T<:Real

    x = similar(ϕ)
    y = similar(λ)
    
    tmp = @. latlon2xy(ϕ, λ; ϕc=ϕc, λ₀=λ₀, 𝑎=𝑎, 𝑒=𝑒)
    foreach(pairs(tmp)) do (i, V)
        x[i] = V[1]
        y[i] = V[2]
    end
    
    return x, y
end
# ----/

"""
Transform x, y data from polar steregraphic system to latitude and longitude

julia> x, y = xy2latlon(x, y, ϕc, λ₀)

Where:
* x : X-axis coordinate of map, scalar or vector [m]
* y : Y-axis coordinate of map, scalar or vector [m]

Optional input parameters:
* ϕ_c: latitude of true scale in degree North, standard parallel (default is -70)
* λ₀: meridian in degree East along positive Y-axis of map (default is 0)
* 𝑎 : Earth radius defined in the projection (default 6378137.0 m, WGS84)
* 𝑒 : eccentricity (default 0.08181919)
OUTPUT:
* ϕ : latitude in degree North
* λ : longitude in degree East

Polar stereographic projection used for satellite polar sea ice studies.
Equations adapted from "Map Projections - A Working Manual! by J.P. Snyder (1987)

"""
function xy2latlon(x::T, y::T; ϕₖ=70, λ₀=0, 𝑎=6378137.0, 𝑒=0.08181919) where T<:Real

    ϕₖ *= π/180
    λ₀ *= π/180
    
    t(α) = tan(π/4-α/2)/√((1-𝑒*sin(α))/(1+𝑒*sin(α)))^𝑒        # Eq. (15-9)  
    m(α) = cos(α)/√(1-(𝑒*sin(α))^2)    # Eq. (14-15)
    tₖ = t(ϕₖ)
    mₖ = m(ϕₖ)
    
    ρ = √(x^2 + y^2)       # Eq. (20-18)
    
    t = ρ*tₖ/mₖ/𝑎          # Eq. (21-40)
    χ = π/2 - 2atan(t)     # Eq. (7-13)

    # Eq. (3-5)
    ϕ = χ + (𝑒^2/2 + 5𝑒^4/24 + 𝑒^6/12 + 13𝑒^8/360)sin(2χ) +
        (7𝑒^4/48 + 29𝑒^6/240 + 811𝑒^8/11520)sin(4χ) +
        (7𝑒^6/120 + 81𝑒^8/1120)sin(6χ) +
        (4279𝑒^8/161280)sin(8χ)
    
    λ = λ₀ + atan(x, -y)        # Eq. (20-16)
    
    # adapting longitude ∈ [-π, π]
    λ = mod(λ+π, 2π) - π
    
    return rad2deg(ϕ), rad2deg(λ)
end

function xy2latlon(x::Array{T}, y::Array{T}; ϕₖ::T=70.0, λ₀::T=0.0, 𝑎=6378137.0, 𝑒=0.08181919) where T<:Real

    
    ϕ = similar(x)
    λ = similar(y)

    tmp = @. xy2latlon(x, y, ϕₖ=ϕₖ, λ₀=λ₀, 𝑎=𝑎, 𝑒=𝑒)
    foreach(pairs(tmp)) do (i, V)
        ϕ[i] = V[1]
        λ[i] = V[2]
    end
    
    return ϕ, λ
end
# ----/


"""
Function to mimic MATLAB meshgrid

julia> X, Y = couple2grid(x, y)

Where:
* x : Vector
* y : Vector
Output:
* X : Matrix (ny, nx) with x repeated along the other dimension,
* Y : Matrix (ny, nx)

where nx = length(x) and ny = length(y)
"""
function couple2grid(x::Vector{T}, y::Vector{T}) where T<:Real
    nx = length(x)
    ny = length(y)

    X = ones(ny)'.*x;
    Y = y'.*ones(nx)
    
    return X, Y
end
# ----/
# end of file
