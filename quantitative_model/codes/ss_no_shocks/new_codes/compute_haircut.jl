include("trandyn_functions.jl")

i_initial = 2
i_final   = 1

B_initial, B_final = B_demand[i_initial], B_demand[i_final]
hairtcut_coeff     = B_final/B_initial

τ_ss     = τvec_comp_an[i_final]
i_τss    = find(nag.τgrid .== τ_ss)[1]
Vf_input = nag.Vf[:,:,:,i_τss]
bp_input = nag.bprime_mat[:,:,:,i_τss]


capT     = 40
max_iter = 5
guess_q  = q_out[i_final] * ones(capT+1)

damp_td_vec = 0.95 * ones(capT)

Λ0 = λast_out[:,i_initial]

Λ0_mat  = reshape(Λ0, length(nag.b_grid_fine), nag.nl)
Λb0     = sum(Λ0_mat,2)
# Λ0_mat_new  = zeros(Λ0_mat)
Λb_new      = zeros(Λb0)
ww_r        = zeros(length(nag.b_grid_fine))

min_ind_vec = Array(Int,length(nag.b_grid_fine))

for i_bf=1:length(nag.b_grid_fine)
    bnew        = nag.b_grid_fine[i_bf]*hairtcut_coeff
    ind_r       = min( max( searchsortedfirst(nag.b_grid_fine, bnew), 2), length(nag.b_grid_fine))
    # ww_r        = ( bnew - nag.b_grid_fine[ind_r-1] ) / ( nag.b_grid_fine[ind_r] - nag.b_grid_fine[ind_r-1] )
    ww_r[i_bf]  = ( bnew - nag.b_grid_fine[ind_r-1] ) / ( nag.b_grid_fine[ind_r] - nag.b_grid_fine[ind_r-1] )

    # min_val, min_ind  = findmin(abs(bnew - nag.b_grid_fine))
    # min_ind_vec[i_bf] = min_ind
    min_ind_vec[i_bf] = ind_r
end

for i_bf=1:length(nag.b_grid_fine)
    # term1 = sum(Λ0_mat[min_ind_vec[min_ind_vec.==i_bf], :],1)*ww_r[i_bf]            # when it is to the right
    # term2 = sum(Λ0_mat[min_ind_vec[(min_ind_vec-1).==i_bf], :],1)*(1-ww_r[i_bf])    # when it is to the left
    # Λ0_mat_new[i_bf,:] = term1 + term2

    # Λ0_mat_new[i_bf,:] = sum(Λ0_mat[min_ind_vec[min_ind_vec.==i_bf], :],1)

    Λb_new[i_bf] = sum(Λb0[min_ind_vec.==i_bf].*ww_r[min_ind_vec.==i_bf]) + sum(Λb0[(min_ind_vec-1).==i_bf].*(1-ww_r[(min_ind_vec-1).==i_bf]))

    # Λb_new[i_bf] = sum(Λb0[min_ind_vec.==i_bf])

end

_, Πl_inv   = eigs(nag.Πl', nev=1)
Πl_inv      = real(squeeze(Πl_inv, 2))
Πl_inv      = Πl_inv ./ sum(Πl_inv)

Λ0_new = kron(Πl_inv, Λb_new)

# fig,ax = subplots()
# ax[:plot](nag.b_grid_fine,[Λb0 Λb_new])





q_dyn_mat_out, q_dyn_out, τ_dyn, Λ_trandyn, bp_cell_temp, qgrid_cell, Vft_cell = solve_trandyn(nag, i_τss, capT, Λ0_new, guess_q, max_iter, damp_td_vec, B_demand[i_final])


bpT = bp_cell_temp[1]

itp_bT = interpolate((nag.bgrid,nag.lgrid,qgrid_cell[1]),bpT,Gridded(Linear()))
itp_bp = interpolate((nag.bgrid,nag.lgrid,nag.qgrid),nag.bprime_mat[:,:,:,i_τss],Gridded(Linear()))

bT_mat = zeros(nag.nb, nag.nl)
bp_mat = zeros(nag.nb, nag.nl)

for (i_b, b_v) in enumerate(nag.bgrid)
    for (i_l,l_v) in enumerate(nag.lgrid)
        bT_mat[i_b, i_l] = itp_bT[b_v,l_v,q_out[i_final]]
        bp_mat[i_b, i_l] = itp_bp[b_v,l_v,q_out[i_final]]
    end
end

maxabs(bT_mat - bp_mat)
λast_out_ee = λast_out[:,i_final]

valbT, valbp = 0.0, 0.0
for i_s=1:size(nag.snodes,1)
    valbT += itp_bT[nag.snodes[i_s,1],nag.snodes[i_s,2],q_out[i_final]]*λast_out_ee[i_s]
    valbp += itp_bp[nag.snodes[i_s,1],nag.snodes[i_s,2],q_out[i_final]]*λast_out_ee[i_s]
end

τ_out = nag.e_bar + (1.0 - q_dyn_out)*B_final - nag.τY



# fig, ax = subplots()
# ax[:plot](nag.b_grid_fine, [sum(reshape(Λ_trandyn[:,end], length(nag.b_grid_fine), nag.nl),2) sum(reshape(λast_out[:,i_final], length(nag.b_grid_fine), nag.nl),2)])
# ax[:set_title]("Convergence of Bonds Distribution")






""" Check Solution """
# itp_Vftcell = interpolate((nag.bgrid, nag.lgrid, qgrid_cell[1]), Vft_cell[1], Gridded(Linear()))
# itp_Vftemp  = interpolate((nag.bgrid, nag.lgrid, nag.qgrid), nag.Vf[:, :, :, i_τss], Gridded(Linear()))

# vfcell_temp = zeros(nag.nb, nag.nl)
# vf_temp     = zeros(nag.nb, nag.nl)

# for (i_b, b_v) in enumerate(nag.bgrid)
#     for (i_l, l_v) in enumerate(nag.lgrid)
#         vfcell_temp[i_b, i_l]   = itp_Vftcell[b_v, l_v, q_out[1]]
#         vf_temp[i_b, i_l]       = itp_Vftemp[b_v, l_v, q_out[1]]
#     end
# end

# fig,ax = subplots()
# for il = 1:nag.nl
#     ax[:plot](nag.bgrid, vfcell_temp[:, il], color="black")
#     ax[:plot](nag.bgrid, vf_temp[:, il], color="grey")
# end





##### viejo ######
# # """ First run ssnag_main_vq.jl """

# i_initial = 1
# i_final   = 1

# B_initial, B_final = B_demand[i_initial], B_demand[i_final]
# hairtcut_coeff     = B_final/B_initial

# τ_ss     = τvec_comp_an[i_final]
# i_τss    = find(nag.τgrid .== τ_ss)[1]
# Vf_input = nag.Vf[:,:,:,i_τss]
# bp_input = nag.bprime_mat[:,:,:,i_τss]


# capT     = 10
# max_iter = 2
# guess_q  = q_out_ee * ones(capT+1)

# damp_td_vec = 1.0 * ones(capT)


# Λ0 = λast_out[:,i_initial]

# # Λ0 = copy(λast_out_ee)

# q_dyn_mat_out, q_dyn_out, τ_dyn, Λ_trandyn, bp_mat_out, Vf_trandyn_out = solve_trandyn(nag, i_final, capT, hairtcut_coeff, Λ0, guess_q, max_iter, damp_td_vec, B_final);


# dens_b_trandyn = Array(Float64,length(nag.b_grid_fine),capT+1)
# dens_b_trandyn[:,1] = sum(reshape(Λ0,length(nag.b_grid_fine),nag.nl),2)
# for tt=2:capT+1
#     dens_b_trandyn[:,tt] = sum(reshape(Λ_trandyn[:,tt],length(nag.b_grid_fine),nag.nl),2)
# end

# fig, ax = subplots()
# ax[:plot](1:1:(capT+1),dens_b_trandyn[[1;50;10],:]')

# Vf_trandyn_itp = interpolate((nag.bgrid,nag.lgrid,nag.qgrid),Vf_trandyn_out,Gridded(Linear()))

# welf_trandyn = 0.0
# for i_s=1:size(nag.snodes,1)
#     welf_trandyn += Λ0[i_s]*Vf_trandyn_itp[nag.snodes[i_s,1],nag.snodes[i_s,2],q_dyn_out[1]]
# end

# c_welf_trandyn = (welf_trandyn*(1.0-nag.γ))^(1.0/(1.0-nag.γ))
# c_welf_start   = (welf_out[end]*(1.0-nag.γ))^(1.0/(1.0-nag.γ))

# change_c = (c_welf_trandyn/c_welf_start-1.0)*100









# function solve_backwards(nag::NAG_ECM,q_dyn_input,τ_dyn_input,capT,i_τend)

#     @printf("Solving Backwards\n")

#   """ Vf_end tiene que ser evaluada en el τ del steady-state """

#     # q_dyn_input contains the guess for the path for q_{t}, starting from q(τ0) and ending on q(τ1) for the last two elements.

#     wage = nag.A_bar

#     # We start from the steady state value function.
#     Vfp_temp = nag.Vf[:,:,:,i_τend]
#     bp_end   = nag.bprime_mat[:,:,:,i_τend]
#     bp_cell  = []

#     for ii=0:capT-1
#         trev = capT - ii

#         Vft_new = Array(Float64, nag.nb, nag.nl, nag.nq)
#         bpt_new = Array(Float64, nag.nb, nag.nl, nag.nq)

#         qt1 = q_dyn_input[trev+1]
#       τt  = τ_dyn_input[trev]

#         if qt1<nag.qgrid[1]||qt1>nag.qgrid[end]
#             itp_temp = extrapolate(interpolate((nag.bgrid, nag.lgrid, nag.qgrid), Vfp_temp, Gridded(Linear())),Linear())
#         else
#             itp_temp = interpolate((nag.bgrid, nag.lgrid, nag.qgrid), Vfp_temp, Gridded(Linear()))
#         end
        
#       Vfp_temp_itp = zeros(nag.nb, nag.nl)

#       # Evaluate RHS of value at tomorrow's guess of q
#       for (i_l,l_v) in enumerate(nag.lgrid)
#           for (i_b, b_v) in enumerate(nag.bgrid)
#               Vfp_temp_itp[i_b, i_l] = itp_temp[b_v, l_v, qt1]
#           end
#       end

#       EVf_trandyn = zeros(nag.nb, nag.nl)
#       compute_EV_trandyn!(nag,Vfp_temp_itp, EVf_trandyn)

#         norm(EVf_trandyn) > 0.0 || throw(error("equal matrices at solve_backwards"))

#       for (i_q, q_v) in enumerate(nag.qgrid)  # Loop over q_t
#           for (i_l, l_v) in enumerate(nag.lgrid)
#               for (i_b, b_v) in enumerate(nag.bgrid)

#                     rhs0 = -1e14
#                     bpol = -1.0

#                     for (i_bp, bp_v) in enumerate(nag.bgrid)

#                         cc    = wage*l_v*(1.0-nag.τY) + b_v  - τt - bp_v*q_v
#                         if cc <=0.0 ; rhs = -1e14
#                         else
#                             rhs   = ufun(nag,cc) + nag.β*EVf_trandyn[i_bp, i_l]
#                         end

#                         if rhs>rhs0
#                             rhs0 = rhs
#                             bpol = bp_v
#                         end

#                     end

#                     Vft_new[i_b, i_l, i_q] = rhs0
#                   bpt_new[i_b, i_l, i_q] = bpol

#               end
#           end
#       end

#         maximum(bpt_new)>0.0 ||throw(error("bp too low"))

#         if ii == 0
#             copy!(Vfp_temp,Vfp_temp_itp)
#             # push!(bp_cell,bp_end)
#             push!(bp_cell, bpt_new)
#         else
#             push!(bp_cell, bpt_new)
#             copy!(Vfp_temp, Vfp_temp_itp)
#         end

        
#     end

#     return bp_cell, Vfp_temp
# end

