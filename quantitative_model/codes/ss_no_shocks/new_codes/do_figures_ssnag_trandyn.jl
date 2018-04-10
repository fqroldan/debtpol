

using PyPlot

fig,ax = subplots()
ax[:plot]([τ_out ones(capT+1)*τvec_comp_an[1] ones(capT+1)*τvec_comp_an[2]])


fig,ax = subplots()
ax[:plot](1./q_dyn_out-1);ax[:set_ylim]([1./nag.qgrid[end]-1,1/nag.qgrid[1]-1])
ax[:plot](ones(q_dyn_mat_out[:,1])./q_out[i_initial] -1, color="grey", label="initial q")
ax[:plot](ones(q_dyn_mat_out[:,1])./q_out[i_final] -1, color="black", label="final q")
ax[:legend](loc=4);
ax[:set_title]("Equilibrium Interest Rate During Transition");ax[:set_xlabel]("Time")

fig,ax = subplots()
iter = 3
for ii=size(q_dyn_mat_out,2)-iter:size(q_dyn_mat_out,2)
    ax[:plot](1./q_dyn_mat_out[:,ii]-1, color="blue", alpha=ii/iter)#;ax[:set_ylim]([1./nag.qgrid[end]-1,1/nag.qgrid[1]-1])
end
ax[:plot](ones(q_dyn_mat_out[:,1])./q_out[i_initial] -1, color="grey", label="initial q")
ax[:plot](ones(q_dyn_mat_out[:,1])./q_out[i_final] -1, color="black", label="final q")
ax[:legend](loc=1);
ax[:set_title]("Interest Rate Convergence");ax[:set_xlabel]("Time")




Vf_trandyn_itp = interpolate((nag.bgrid,nag.lgrid,qgrid_cell[end]),Vft_cell[end],Gridded(Linear()))

welf_trandyn = 0.0
for i_s=1:size(nag.snodes,1)
    welf_trandyn += Λ0_new[i_s]*Vf_trandyn_itp[nag.snodes[i_s,1]*hairtcut_coeff,nag.snodes[i_s,2],q_dyn_out[1]]
end


dens_b_trandyn      = Array(Float64,length(nag.b_grid_fine),capT+1)
dens_b_trandyn[:,1] = sum(reshape(Λ0_new,length(nag.b_grid_fine),nag.nl),2)


val_quantiles     = Array(Float64,capT+1, 10)
i_qs = [collect(1:9); 9.9]
for (iquant, quant) in enumerate(i_qs)
    val_quantiles[1,iquant]  = nag.b_grid_fine[indmin(abs(cumsum(dens_b_trandyn[:,1])-quant/10))]
end

for tt=2:capT+1
    dens_b_trandyn[:,tt] = sum(reshape(Λ_trandyn[:,tt],length(nag.b_grid_fine),nag.nl),2)
    for (iquant, quant) in enumerate(i_qs)
        val_quantiles[tt,iquant] = nag.b_grid_fine[indmin(abs(cumsum(dens_b_trandyn[:,tt])-quant/10))]
    end
end

fix,ax = subplots()
for (iquant, quant) in enumerate(i_qs)
    ax[:plot](val_quantiles[:,iquant])
end
ax[:set_title]("Quantiles of Bond Holdings");ax[:set_xlabel]("Time")


# Experiment 2: Haircut and Dynamics

welf_hh = Array(Float64,capT+1,3)

i_τss_initial = find(nag.τgrid .== τvec_comp_an[i_initial])[1]

itp_Vf = interpolate((nag.bgrid,nag.lgrid,nag.qgrid),nag.Vf[:,:,:,i_τss_initial],Gridded(Linear()))

i_l_welf = 2#nag.nl
welf_hh[1,1] = itp_Vf[b_poster[1],nag.lgrid[i_l_welf],q_out[i_initial]]
welf_hh[1,2] = itp_Vf[b_poster[2],nag.lgrid[i_l_welf],q_out[i_initial]]
welf_hh[1,3] = itp_Vf[b_poster[3],nag.lgrid[i_l_welf],q_out[i_initial]]

soc_welf = Array(Float64,capT+1)

soc_welf[1] = welf_out[i_initial]


for ii=0:capT-1
    trev = capT - ii

    Vf_temp    = Vft_cell[trev]
    qgrid_ext  = qgrid_cell[trev]
    itp_Vftemp = interpolate((nag.bgrid,nag.lgrid,qgrid_ext),Vf_temp,Gridded(Linear()))

    if q_dyn_out[ii+1]<qgrid_ext[1] || q_dyn_out[ii+1]>qgrid_ext[end]
        itp_Vf_obj = extrapolate(itp_Vftemp,Linear())
    else
        itp_Vf_obj = itp_Vftemp
    end

    welf_hh[ii+2,1] = itp_Vf_obj[b_poster[1],nag.lgrid[i_l_welf],q_dyn_out[ii+1]]#cons_eq(nag,itp_Vf_obj[b_poster[1],1.0,q_dyn_out[ii+1]])
    welf_hh[ii+2,2] = itp_Vf_obj[b_poster[2],nag.lgrid[i_l_welf],q_dyn_out[ii+1]]#cons_eq(nag,itp_Vf_obj[b_poster[2],1.0,q_dyn_out[ii+1]])
    welf_hh[ii+2,3] = itp_Vf_obj[b_poster[3],nag.lgrid[i_l_welf],q_dyn_out[ii+1]]#cons_eq(nag,itp_Vf_obj[b_poster[3],1.0,q_dyn_out[ii+1]])

    welf_soc_val = 0.0
    for i_s=1:nag.num_nodes
        welf_soc_val += itp_Vf_obj[nag.snodes[i_s,1],nag.snodes[i_s,2],q_dyn_out[ii+1]] * Λ_trandyn[i_s, ii+1]
    end

    soc_welf[ii+2] = welf_soc_val
end

fig,ax = subplots()
for ii=1:3
    ax[:plot](welf_hh[:,ii], label="b=$(b_poster[ii])");ax[:legend](loc=4);ax[:set_title]("Implications of a Haircut on Welfare")
end
ax[:set_xlabel]("Time")

fig,ax = subplots()
ax[:plot](soc_welf)
;ax[:set_title]("Aggregate Welfare Upon Haircut")
