function ∫fdh(x::AbstractArray, H::AbstractVector; h₀=Real[], hₜ=Real[])
    # getting dimensions:
    nheight, ntime = size(x)

    # checking dimensions:
    nheight != length(H) && @error "Heights muss have same length as matrix 1st dim."
    !isempty(h₀) && ntime != length(h₀) && @error "Optional variable h₀ must have length same as matrix 2nd dim."
    !isempty(hₜ) && ntime != length(hₜ) && @error "Optional variable hₜ must have length same as matrix 2nd dim."
	
    # defining the derivate of H
    δh = Vector{eltype(x)}(undef, length(H)) .= 0
    δh[2:end] = diff(H)
    δh[1] =  δh[2:end] |> minimum

    # integrating variable within limits:
    𝐼₀ₜ = fill(NaN32, ntime)
    for (i, X) ∈ enumerate(eachcol(x))
        # finding the limits for integration:
        i0 = isempty(h₀) ? 1 : !isnan(h₀[i]) ? argmin(abs.(H .- h₀[i])) : nothing
        it = isempty(hₜ) ? nheight : !isnan(hₜ[i]) ? argmin(abs.(H .- hₜ[i])) : nothing

        isnothing(i0) && continue
        isnothing(it) && continue

	lims = i0:it

        Xh = let Xh=X[lims] 
            tmp = ismissing.(Xh) .|| isnan.(Xh)
            all(tmp) && continue
            Xh[tmp] .= 0
            Xh
        end
        #Xh = (collect∘skipmissing)(Xh)
        
	𝐼₀ₜ[i] = Xh'*vec(δh[lims]) 

        
    end
    return 𝐼₀ₜ
end
#----/

