function fill_sparse_nag(grid::Vector, pol::Vector)

    # This function takes as input the vector b_{i} and b_prime_{i}
    # and yields the needed values for constructing the sparse matrix.

    #== Initialize vectors ==#
    I_vec   = collect( 1:length(pol) )
    J_vec   = similar(I_vec)

    w_right = Array{Real}(length(pol))


    for (i,value) in enumerate(pol)
        indj   = min(max( searchsortedfirst(grid, value), 2),length(grid))  ## NOTE: 1st value ≥ in grid

        J_vec[i]   = indj
        value      = max(min(value,grid[end]),grid[1])
        w_right[i] = ( value - grid[indj-1] ) / ( grid[indj] - grid[indj-1] )
    end

    #== Repeat indices ==#
    I_vec   = [I_vec; I_vec ]
    J_vec   = [J_vec; J_vec-1 ]
    weights = [w_right; 1.0 - w_right]

    return I_vec, J_vec, weights
end

function get_sparse_nag(nag::NAG_ECM, qval, qgrid_input, bpmat_input)
    # In this function we generate a sparse transition matrix Π_{t},
    # such that Λ_{t+1}=Π_{t}'Λ_{t}.
    # Π_{t} = sparse(I,J,V), where I is a vector indicating the location of the row,
    # J is a vector containing the location of the column, and V is a vector with 
    # the corresponding values to fill in.
    # We will actually be computing Π_{t}'.

    ln_mk    = length(nag.lgrid)
    nb_fine  = length(nag.b_grid_fine)

    #== Initialize Vectors for sparse matrix ==#
    Ic      = Array{Int64}(0)
    Ir      = Array{Int64}(0)
    V_vec   = Array{Real}(0)         # should contain Reals

    itp_obj_bi = interpolate((nag.bgrid, nag.lgrid, qgrid_input), bpmat_input, Gridded(Linear()))

    for i_l = 1:ln_mk
        # Generate matrix of next period's bonds policy given current state
        bprime = Array{Float64}(nb_fine)
        for i_b=1:nb_fine
            bprime[i_b] = itp_obj_bi[nag.b_grid_fine[i_b], nag.lgrid[i_l], qval]
        end

        # Here we change the output to the transpose
        ic, ir, vv = fill_sparse_nag(collect(nag.b_grid_fine), bprime)

        for i_l1 = 1:ln_mk
            prob = nag.Πl[i_l, i_l1] # Get Markov transition probabilities
            ioff = (i_l-1)*nb_fine   # Generate index
            joff = (i_l1-1)*nb_fine  # Generate index

            #== Append on indices ==#
            append!(Ic, ic + ioff)
            append!(Ir, ir + joff)
            append!(V_vec, prob * vv )
        end
    end
    nΠaggr = nb_fine*ln_mk

    return sparse(Ir, Ic, V_vec, nΠaggr, nΠaggr)
end


function compute_SS_once(nag::NAG_ECM, qval, qgrid_input, bpmat_input)

    # Compute aggregate demand for bonds and invariant distribution at a given price
    
    NS      = size(nag.snodes,1)

    itp_obj_bi = interpolate((nag.bgrid, nag.lgrid, qgrid_input), bpmat_input, Gridded(Linear()))

    # Comute the transition operator for the histogram
    Πst       = get_sparse_nag(nag, qval, qgrid_input, bpmat_input)
    _, eigenv = eigs(Πst, nev = 1);
    eigenv    = real( squeeze(eigenv,2) )
    λast      = eigenv ./ sum(eigenv);

    #== check all elements are positive ==#
    minimum(λast) > - 1e-8 || throw(error("Negative histogram weight %.8f",minimum(λast)))
    λast[λast .< 1e-12] = 0.0
    λast = λast ./ sum(λast)

    #== CHECK consistency ==#
    err_invDist = maximum( abs.( Πst*λast - λast) );
    # @printf("   - Invariant distribution Converged to %.2e \n", err_invDist)
    err_invDist<1e-8 || throw(error("Invariant distribution is incorrect"))


    # Compute the aggregates
    B_demand_out = 0.0

    for i_s=1:NS
        bprime        = itp_obj_bi[nag.snodes[i_s,1], nag.snodes[i_s,2], qval]
        B_demand_out += bprime*λast[i_s]
    end

    return λast, B_demand_out[1]
end


function compute_bond_supply(nag::NAG_ECM, τval)

    # Compute aggregate supply of bonds from GG b.c., for each possible q_v

    e_ss = e_fun(nag,nag.A_bar)
    wage = nag.A_bar
    
    tax_coll = 1.0 - τval * dot(wage*nag.lgrid.^(1.0-nag.τ2), nag.Πl_ast)

    Bs_mat_out = Array{Float64}(nag.nq)

    for (i_q,q_val) in enumerate(nag.qgrid)

        Bs_mat_out[i_q] = (e_ss - nag.τY*wage*nag.L_bar - τval)/(q_val - 1.0)
        if nag.prog_tax
            Bs_mat_out[i_q] = (e_ss - tax_coll)/(q_val - 1.0)
        end
        Bs_mat_out[i_q] >= 0.0 || throw(error("Negative B at compute_bond_supply"))
    end

    return Bs_mat_out
end

function compute_ss_root_once(nag::NAG_ECM, τval, qgrid_input, bpmat_input)

    # Find q such that AD = AS of bonds

    e_ss = e_fun(nag, nag.A_bar)
    wage = nag.A_bar
    tax_coll = 1.0 - τval * dot(wage*nag.lgrid.^(1.0-nag.τ2), nag.Πl_ast)

    λast_out      = Array{Float64}(size(nag.snodes,1))
    B_demand_temp = 0.0
    B_supply_temp = 0.0

    function temp_fun(q_input)
        λast_out, B_demand_temp = compute_SS_once(nag, q_input, qgrid_input, bpmat_input)

        B_supply_temp = (e_ss - nag.τY*wage*nag.L_bar - τval)/(q_input - 1.0)
        if nag.prog_tax
            B_supply_temp = (e_ss - tax_coll)/(q_input - 1.0)
        end
        B_supply_temp >= 0.0 || throw(error("Negative B at compute_ss_root"))

        dif_B = abs(B_supply_temp - B_demand_temp)
        
        return dif_B
    end

    res   = optimize(temp_fun, nag.qgrid[1], nag.qgrid[end], GoldenSection(), iterations = 100) # 500
    q_out = Optim.minimizer(res)[1]

    return λast_out, B_demand_temp, B_supply_temp, q_out
end

function compute_EV_ss!(nag::NAG_ECM, Vfp_input, EVf_ss)

    for i_l=1:nag.nl, i_bp=1:nag.nb
        EV = 0.0
        for i_lp=1:nag.nl
            EV += nag.Πl[i_l, i_lp]*Vfp_input[i_bp, i_lp]
        end
        EVf_ss[i_bp, i_l] = EV
    end

    return Void
end

function iterate_VF_ss!(nag::NAG_ECM, Vf_input, bpol_input, cpol_input, qgrid_ext, itp_EVf_input, q_eqm)

    iter = 0
    dist = 10.0

    VF_old = zeros(nag.nb, nag.nl, length(qgrid_ext))

    while dist > 1.0e-5 && iter < 500

        iter += 1

        for (i_q, q_v) in enumerate(qgrid_ext)  # Loop over q_t on extended qgrid.
            for (i_l, l_v) in enumerate(nag.lgrid)
                for (i_b, b_v) in enumerate(nag.bgrid)
                    Vf_input[i_b, i_l, i_q] = ufun(nag, cpol_input[i_b, i_l, i_q]) + nag.β*itp_EVf_input[bpol_input[i_b, i_l, i_q], l_v]
                end
            end
        end

        ### Update expetation of VF

        itp_Vfp_temp = interpolate((nag.bgrid, nag.lgrid, qgrid_ext), Vf_input, Gridded(Linear()))

        # if q_eqm<qgrid_ext[1] || q_eqm>qgrid_ext[end]
        #     itp_Vfp_obj = extrapolate( itp_Vfp_temp, Linear() )
        # else
            # itp_Vfp_obj = itp_Vfp_temp
        # end

        Vfp_temp_itp = Array{Float64}( nag.nb, nag.nl)

        # Evaluate RHS of value at tomorrow's guess of q
        for (i_l,l_v) in enumerate(nag.lgrid)
            for (i_b, b_v) in enumerate(nag.bgrid)
                # Vfp_temp_itp[i_b, i_l] = itp_Vfp_obj[b_v, l_v, q_eqm]
                Vfp_temp_itp[i_b, i_l] = itp_Vfp_temp[b_v, l_v, q_eqm]
            end
        end    

        EVf_ss = zeros(nag.nb, nag.nl)
        compute_EV_ss!(nag, Vfp_temp_itp, EVf_ss)

        itp_EVf_input = interpolate((nag.bgrid, nag.lgrid), EVf_ss, Gridded(Linear())) 

        dist   = maximum(abs, VF_old - Vf_input)
        copy!(VF_old, Vf_input)
    end

    @printf("Finished SS VF iteration with distance %.5f and %d iterations\n", dist, iter)

    return Void
end

function compute_Vb_once_ss!(nag::NAG_ECM, Vbp_red_input, qgrid_ext, q_eqm, τval, bpt_new, Vbt_new, cpol_new)


    # Compute the expectations
    EVb_ss = zeros(nag.nb, nag.nl)
    compute_EV_ss!(nag, Vbp_red_input, EVb_ss)

    # Create interpolations for EVb, evaluated at q_{t+1}        
    itp_EVbp = interpolate((nag.bgrid, nag.lgrid), EVb_ss, Gridded(Linear()))

    # Compute new policies and new envelope        
    for (i_q, q_v) in enumerate(qgrid_ext)  # Loop over q_t on extended qgrid.
        for (i_l, l_v) in enumerate(nag.lgrid)
            τl = 1.0 - τ_v*l_v^(-nag.τ2)
            for (i_b, b_v) in enumerate(nag.bgrid)

                # First check whether the constraint binds
                bp_v = nag.bgrid[1]
                euler_val = uprime(nag, nag.A_bar*l_v*(1.0-nag.τY) + b_v - q_v*bp_v - τval) * q_v - nag.β * itp_EVbp[bp_v, l_v]
                if nag.prog_tax
                    euler_val = uprime(nag, nag.A_bar*l_v*(1.0-τl) + b_v - q_v*bp_v - τval) * q_v - nag.β * itp_EVbp[bp_v, l_v]
                end
                if euler_val >= 0.0
                    bpol = nag.bgrid[1]
                else
                    # If constraint not binding, find b' that solves the Euler equation at equality
                    l_bound = nag.bgrid[1]
                    u_bound = max((nag.A_bar*l_v*(1.0-nag.τY) + b_v  - τval)/q_v - 1e-5, nag.bgrid[1]+1e-5)
                    if nag.prog_tax
                        u_bound = max((nag.A_bar*l_v*(1.0-τl) + b_v  - τval)/q_v - 1e-5, nag.bgrid[1]+1e-5)
                    end
                    if nag.prog_tax
                        res     = optimize((bp_v -> abs(uprime(nag,nag.A_bar*l_v*(1.0-τl) + b_v - q_v*bp_v - τval)*q_v - nag.β*itp_EVbp[bp_v,l_v])), l_bound, u_bound, GoldenSection(), iterations = 15)
                    else
                        res     = optimize((bp_v -> abs(uprime(nag,nag.A_bar*l_v*(1.0-nag.τY) + b_v - q_v*bp_v - τval)*q_v - nag.β*itp_EVbp[bp_v,l_v])), l_bound, u_bound, GoldenSection(), iterations = 15)
                    end
                    bpol    = Optim.minimizer(res)[1]
                    # bpol    = fzeros((bp_v -> uprime(nag,nag.A_bar*l_v*(1.0-nag.τY) + b_v - q_v*bp_v - τval)*q_v - nag.β*itp_EVbp[bp_v,l_v]), l_bound, u_bound)[1]
                end

                # Deduce consumption from budget constraint
                if nag.prog_tax
                    cpol = max(nag.A_bar*l_v*(1.0-τl) + b_v - q_v*bpol - τval, 1e-14)
                else
                    cpol = max(nag.A_bar*l_v*(1.0-nag.τY) + b_v - q_v*bpol - τval, 1e-14)
                end
                cpol>0.0 || throw(error("negative consumption $(round(cpol,3)) in solve_backwards_ECM with bonds $(round(bpol,3)) and Euler $(round(euler_val,3))"))

                # Save current choices and value functions
                bpt_new[i_b, i_l, i_q]  = bpol
                Vbt_new[i_b, i_l, i_q]  = uprime(nag, cpol)  # Exploit envelope
                cpol_new[i_b, i_l, i_q] = cpol
            end
        end
    end

    return Void

end

function reduce_dim_mat_ss(nag::NAG_ECM, mat_input, qgrid_ext, q_eqm)

    mat_temp_itp = Array{Float64}( nag.nb, nag.nl)  

    # Interpolate to compute old Vb and old policy
    itp_mat_temp = interpolate((nag.bgrid, nag.lgrid, qgrid_ext), mat_input, Gridded(Linear()))

    # if q_eqm<qgrid_ext[1] || q_eqm>qgrid_ext[end]
    #     itp_mat_obj = extrapolate( itp_mat_temp, Linear() )
    # else
        # itp_mat_obj = itp_mat_temp
    # end

    # Evaluate RHS of value at tomorrow's guess of q
    for (i_l,l_v) in enumerate(nag.lgrid)
        for (i_b, b_v) in enumerate(nag.bgrid)
            # mat_temp_itp[i_b, i_l] = itp_mat_obj[b_v, l_v, q_eqm]
            mat_temp_itp[i_b, i_l] = itp_mat_temp[b_v, l_v, q_eqm]
        end
    end

    return mat_temp_itp

end


function iterate_Vb_ss(nag::NAG_ECM, Vb_red_input, bp_red_input, qgrid_ext, q_eqm, τval)

    iter = 0
    dist = 10.0
    dist_Vb, dist_bp = 10.0, 10.0


    Vb_input    = Array{Float64}( nag.nb, nag.nl, length(qgrid_ext) )
    bpmat_input = Array{Float64}( nag.nb, nag.nl, length(qgrid_ext) )
    cpol_new    = Array{Float64}( nag.nb, nag.nl, length(qgrid_ext) )

    while dist > 1e-5 && iter < 300

        iter += 1

        Vbt_new  = Array{Float64}( nag.nb, nag.nl, length(qgrid_ext) )
        bpt_new  = Array{Float64}( nag.nb, nag.nl, length(qgrid_ext) )        

        compute_Vb_once_ss!(nag, Vb_red_input, qgrid_ext, q_eqm, τval, bpt_new, Vbt_new, cpol_new)

        if iter == 1

            Vb_red_new = reduce_dim_mat_ss(nag, Vbt_new, qgrid_ext, q_eqm)
            bp_red_new = reduce_dim_mat_ss(nag, bpt_new, qgrid_ext, q_eqm)

            dist_Vb = maximum(abs, Vb_red_new - Vb_red_input)
            dist_bp = maximum(abs, bp_red_new - bp_red_input)

        else

            dist_Vb = maximum(abs, Vbt_new - Vb_input)
            dist_bp = maximum(abs, bpt_new - bpmat_input)
        end

        dist = max(dist_Vb, dist_bp)

        Vb_input    = copy(Vbt_new)
        bpmat_input = copy(bpt_new)

        Vb_red_input = reduce_dim_mat_ss(nag, Vb_input,    qgrid_ext, q_eqm)
        bp_red_input = reduce_dim_mat_ss(nag, bpmat_input, qgrid_ext, q_eqm)
        
    end

    @printf("Vb # iterations %d, dist_Vb %.3g, dist_bp %.3g\n", iter, dist_Vb, dist_bp)

    return (Vb_input, bpmat_input, cpol_new)
end


function compute_ss_root(nag::NAG_ECM, τval; do_simple = true)

    iter_ss      = 0
    q_old        = 0.0
    B_old        = 0.0
    dist_r_ss    = 10.0
    dist_B_ss    = 10.0

    i_τ = find(nag.τgrid .== τval)[1]

    λast_out, B_demand_temp, B_supply_temp, q_out = [], [], [], []

    qgrid_ext   = copy(nag.qgrid)
    bpmat_input = copy(nag.bprime_mat[:, :, :, i_τ])
    Vft_input   = copy(nag.Vf[:, :, :, i_τ])
    Vb_input    = copy(nag.Vb_mat[:, :, :, i_τ])
    cpol_mat    = copy(nag.c_mat[:, :, :, i_τ])

    if do_simple
        λast_out, B_demand_temp, B_supply_temp, q_out = compute_ss_root_once(nag, τval, qgrid_ext, bpmat_input)
    else
        while (dist_r_ss > 1e-5 || dist_B_ss>1e-5) && iter_ss<10

            iter_ss += 1

            λast_out, B_demand_temp, B_supply_temp, q_out = compute_ss_root_once(nag, τval, qgrid_ext, bpmat_input)

            # Evaluate Vb and gb at original grid
            Vb_red_initial = reduce_dim_mat_ss(nag, Vb_input,    qgrid_ext, q_out)
            bp_red_initial = reduce_dim_mat_ss(nag, bpmat_input, qgrid_ext, q_out)

            # Create new grid for q
            qgrid_ext = zeros(length(nag.qgrid)+1)
            if q_out > maximum(nag.qgrid)
                qgrid_ext[1:end-1]  = nag.qgrid[:]
                qgrid_ext[end]      = q_out
            elseif q_out < minimum(nag.qgrid)
                qgrid_ext[2:end]    = nag.qgrid[:]
                qgrid_ext[1]        = q_out
            else
                ind_qt1_r = searchsortedfirst(nag.qgrid, q_out)
                ind_qt1_l = ind_qt1_r-1

                if q_out != nag.qgrid[ind_qt1_r] # in this case we have to add the point

                    qgrid_ext[1:ind_qt1_l]     = nag.qgrid[1:ind_qt1_l]
                    qgrid_ext[ind_qt1_l+1]     = q_out
                    qgrid_ext[ind_qt1_l+2:end] = nag.qgrid[ind_qt1_r:end]
                else
                    qgrid_ext = copy(nag.qgrid)
                end
            end

            # Iterate over Vb
            out_iter_Vb  = iterate_Vb_ss(nag, Vb_red_initial, bp_red_initial, qgrid_ext, q_out, τval)

            Vb_input     = out_iter_Vb[1]
            bpmat_input  = out_iter_Vb[2]
            cpol_mat     = out_iter_Vb[3]

            dist_B_ss = abs(B_demand_temp - B_old)
            dist_r_ss = abs(1.0/q_out - 1.0/q_old)
            B_old     = B_demand_temp        
            q_old     = q_out

            @printf("B = %.7f, r = %.7f, dist_B = %.5f, dist_r = %.5f\n", B_demand_temp, 1.0/q_out, dist_B_ss, dist_r_ss)

            size(bpmat_input, 3) == length(qgrid_ext) || throw(error("mirar"))
        end
    end

    # Use new policies to compute new value function (by iteration)

    Vf_red_input = reduce_dim_mat_ss(nag, Vft_input, nag.qgrid, q_out)
    EVf_ss       = zeros(nag.nb, nag.nl)
    compute_EV_ss!(nag, Vf_red_input, EVf_ss)

    itp_EVfp = interpolate((nag.bgrid, nag.lgrid), EVf_ss, Gridded(Linear())) 
    Vft_new  = Array{Float64}( nag.nb, nag.nl, length(qgrid_ext))
    iterate_VF_ss!(nag, Vft_new, bpmat_input, cpol_mat, qgrid_ext, itp_EVfp, q_out)


    if do_simple == false
        λast_out, B_demand_temp, B_supply_temp, q_out = compute_ss_root_once(nag, τval, qgrid_ext, bpmat_input)
        @printf("Finished iteration on SS with %d iterations and dist_1 %.3f and dist_B %.3f\n", iter_ss, dist_r_ss, dist_B_ss)
    end

    out_objects_SS = (λast_out, B_demand_temp, B_supply_temp, q_out, qgrid_ext, Vft_new, bpmat_input, Vb_input)

    return out_objects_SS
end




function compute_welfare_nag(nag::NAG_ECM, Λ, qval, qgrid_input, Vf_input)

    itp_Vf = interpolate((nag.bgrid, nag.lgrid, qgrid_input), Vf_input, Gridded(Linear()))

    NS = size(Λ,1)
    welf_out = 0.0
    for i_s=1:NS
        welf_out += itp_Vf[nag.snodes[i_s,1], nag.snodes[i_s,2], qval] * Λ[i_s]
    end
    return welf_out
end