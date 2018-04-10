using Distributions
function tauchen_fun(ny::Int64, ρ::Float64, σe::Float64; m=3, mu=0.0)
    σy = σe/sqrt((1-ρ^2))
    λ1 = -m*σy; λn = m*σy
    λgrid = linspace(λ1,λn,ny)
    ww = λgrid[2] - λgrid[1]

    distrib = Normal(0,1)

    Π = Array{Float64}(ny,ny)

    for ii=1:ny
        Π[ii,1] = cdf(distrib,(λ1+ww/2-(1-ρ)*mu-ρ*λgrid[ii])/σe) # For j=1
        Π[ii,end] = 1.0-cdf(distrib,(λn-ww/2-(1-ρ)*mu-ρ*λgrid[ii])/σe) # For j=ny
        for jj=2:ny-1
            Π[ii,jj] = cdf(distrib,(λgrid[jj]+ww/2-(1-ρ)*mu-ρ*λgrid[ii])/σe) - cdf(distrib,(λgrid[jj]-ww/2-(1-ρ)*mu-ρ*λgrid[ii])/σe)
        end
    end

    return λgrid,Π
end