
# We are going to solve for the steady state of an economy with idiosyncratic but no aggregate shocks.
# The government can choose default and taxes. Bond prices will come from the evolution of government's debt.

function compute_ctilde_nag!(nag::NAG_ECM, ctilde_mat)
    # This function computes consumption only for the states without default.

    for i_τ=1:nag.nτ,i_q=1:nag.nq,i_l=1:nag.nl,i_b=1:nag.nb
        ctilde_mat[i_b,i_l,i_q,i_τ] = uprime_inv(nag, nag.Vb_mat[i_b,i_l,i_q,i_τ])
    end

    return Void
end

function compute_consumption_nag!(nag::NAG_ECM)

    wage = nag.A_bar

    for (i_τ,τ_v) in enumerate(nag.τgrid)
        for (i_q,q_v) in enumerate(nag.qgrid)
            for (i_l,l_v) in enumerate(nag.lgrid)
                τl = 1.0 - τ_v*l_v^(-nag.τ2)
                for (i_b,b_v) in enumerate(nag.bgrid)
                    nag.c_mat[i_b,i_l,i_q,i_τ]  = wage*l_v*(1.0-nag.τY) + b_v*nag.repay  - τ_v - nag.bprime_mat[i_b,i_l,i_q,i_τ]*q_v
                    if nag.prog_tax
                        nag.c_mat[i_b,i_l,i_q,i_τ]  = wage*l_v*(1.0-τl) + b_v*nag.repay - nag.bprime_mat[i_b,i_l,i_q,i_τ]*q_v
                    end
                end
            end
        end
    end
end

function compute_EVb_nag(nag::NAG_ECM,i_l,bp,qp_v,τp_v,itp_Vb)

    EVb = 0.0

    # if τp_v<nag.τgrid[1]||qp_v<nag.qgrid[1]||bp<nag.bgrid[1]||τp_v>nag.τgrid[end]||qp_v>nag.qgrid[end]||bp>nag.bgrid[end]
    #     itp_obj_Vb    = extrapolate(itp_Vb,Linear())
    # else
        # itp_obj_Vb    = itp_Vb
    # end

    for (i_lp,lp_v) in enumerate(nag.lgrid)
        # EVb += nag.Πl[i_l,i_lp]*itp_obj_Vb[bp,lp_v,qp_v,τp_v]
        EVb += nag.Πl[i_l,i_lp]*itp_Vb[bp,lp_v,qp_v,τp_v]
    end

    return EVb

end

function compute_EV_nag(nag::NAG_ECM,i_l,bp,qp_v,τp_v,itp_Vf)

    EV = 0.0

    # if τp_v<nag.τgrid[1]||qp_v<nag.qgrid[1]||bp<nag.bgrid[1]||τp_v>nag.τgrid[end]||qp_v>nag.qgrid[end]||bp>nag.bgrid[end]
    #     itp_obj_Vf = extrapolate(itp_Vf,Linear())
    # else
        # itp_obj_Vf = itp_Vf
    # end

    for (i_lp,lp_v) in enumerate(nag.lgrid)
        # EV += nag.Πl[i_l,i_lp]*itp_obj_Vf[bp,lp_v,qp_v,τp_v]
        EV += nag.Πl[i_l,i_lp]*itp_Vf[bp,lp_v,qp_v,τp_v]
    end

    return EV

end

function compute_EVf!(nag::NAG_ECM,EVf)
    for i_τ = 1:nag.nτ
        i_τ1 = i_τ
        for i_q = 1:nag.nq
            i_q1 = i_q
            for i_l = 1:nag.nl
                for (i_bp,bp_v) in enumerate(nag.bgrid)

                    EV = 0.0
                    for (i_lp,lp_v) in enumerate(nag.lgrid)
                        EV += nag.Πl[i_l,i_lp]*nag.Vf[i_bp,i_lp,i_q1,i_τ1]
                    end
                    EVf[i_bp,i_l,i_q,i_τ] = EV
                end
            end
        end
    end
    return Void
end

function compute_envelope_nag(nag::NAG_ECM, ctilde_mat)

    # This function computes the bonds policy and the updated guess for the envelope.
    # It takes as inputs the consumption matrices that were previously computed.
    # It also has as input pre-allocated bond matrices that will be overwritten.

    itp_Vb     = interpolate((nag.bgrid,nag.lgrid,nag.qgrid,nag.τgrid),nag.Vb_mat,Gridded(Linear()))
    Vb_mat_new = copy(nag.Vb_mat)

    for (i_τ, τ_v) in enumerate(nag.τgrid)
        τp_v = τ_v

        for (i_q, q_v) in enumerate(nag.qgrid)

            qp_v = q_v

            wage = nag.A_bar

            for (i_l,l_v) in enumerate(nag.lgrid)
                τl = 1.0 - τ_v*l_v^(-nag.τ2)
                for (i_b,b_v) in enumerate(nag.bgrid)

                    bprime = (1.0/q_v)*(wage*l_v*(1.0-nag.τY) + b_v - ctilde_mat[i_b,i_l,i_q,i_τ] - τ_v)
                    if nag.prog_tax
                        bprime = (1.0/q_v)*(wage*l_v*(1.0-τl) + b_v - ctilde_mat[i_b,i_l,i_q,i_τ])
                    end

                    if bprime<0.0
                        nag.bprime_mat[i_b,i_l,i_q,i_τ] = nag.bgrid[1]
                        cc  = wage*l_v*(1.0-nag.τY)+b_v  - τ_v - qp_v*nag.bgrid[1] # cambio aca!! antes tenia 0
                        if nag.prog_tax
                            cc  = wage*l_v*(1.0-τl)+b_v - qp_v*nag.bgrid[1] # cambio aca!! antes tenia 0
                        end
                        Vb_mat_new[i_b,i_l,i_q,i_τ] = uprime(nag,cc)
                    else
                        nag.bprime_mat[i_b, i_l, i_q, i_τ] = bprime
                        EVb = compute_EVb_nag(nag, i_l, bprime, qp_v, τp_v, itp_Vb)
                        Vb_mat_new[i_b, i_l, i_q, i_τ]     = (nag.β / q_v) * EVb
                    end
                end
            end

        end
    end

    return Vb_mat_new
end


function iterate_envelope_nag!(nag::NAG_ECM, ctilde_mat)

    iter = 0
    dist = 10.0
    damp = 0.9


    compute_ctilde_nag!(nag, ctilde_mat)

    Vb_mat_old = copy(nag.Vb_mat)

    while iter<nag.maxit && dist>nag.tol_dist

        iter += 1

        bprime_old = copy(nag.bprime_mat)

        Vb_mat_new = compute_envelope_nag(nag, ctilde_mat)
        compute_ctilde_nag!(nag, ctilde_mat)

        nag.Vb_mat = damp*nag.Vb_mat + (1.0-damp)*Vb_mat_new

        dist  = maximum(abs,(Vb_mat_new[Vb_mat_new.>0.0]-Vb_mat_old[Vb_mat_new.>0.0]) ./ Vb_mat_old[Vb_mat_new.>0.0])
        distb = maximum(abs,bprime_old-nag.bprime_mat)

        # copy!(Vb_mat_old, Vb_mat_new)
        copy!(Vb_mat_old, nag.Vb_mat)

        if iter%50 == 0
            @printf("Iteration # %d with distance in Vb of %.5f and on bonds of %.5f\n",iter,dist,distb)
        end
    end

    return Void
end

function compute_Vf_once_nag!(nag::NAG_ECM,EVf)

    Vf_new = Array{Float64}(nag.nb,nag.nl,nag.nq)

    for i_τ =1:nag.nτ
        for i_q =1:nag.nq
            for i_l=1:nag.nl
                for i_b=1:nag.nb
                    nag.Vf[i_b,i_l,i_q,i_τ] = ufun(nag,nag.c_mat[i_b,i_l,i_q,i_τ]) + nag.β*EVf[i_b,i_l,i_q,i_τ]
                end
            end
        end
    end
end

function compute_Vf_nag!(nag::NAG_ECM)

    iter_vf = 0
    dist_vf = 10.0

    Vf0 = copy(nag.Vf)

    while iter_vf<nag.maxit && dist_vf>nag.tol_dist

        iter_vf += 1

        itp_Vf = interpolate((nag.bgrid,nag.lgrid,nag.qgrid,nag.τgrid),nag.Vf,Gridded(Linear()))

        # Compute the expectation of the value function given state today
        EVf = Array{Float64}(size(nag.Vf))

        for (i_τ, τ_v) in enumerate(nag.τgrid)
            τp_v = τ_v
            for (i_q, q_v) in enumerate(nag.qgrid)
                qp_v = q_v
                for i_l=1:nag.nl
                    for (i_b,b_v) in enumerate(nag.bgrid)
                        bp_v = nag.bprime_mat[i_b,i_l,i_q,i_τ]
                        EVf[i_b,i_l,i_q,i_τ] = compute_EV_nag(nag,i_l,bp_v,qp_v,τp_v,itp_Vf)
                    end
                end
            end
        end

        compute_Vf_once_nag!(nag,EVf)

        dist_vf = maximum(abs,(nag.Vf-Vf0)./Vf0)

        Vf0 = copy(nag.Vf)

        if iter_vf%50 == 0
            @printf("Distance in Vf is %.5f at iteration %d\n",dist_vf,iter_vf)
        end

    end

end


function compute_bonds_in_grid!(nag::NAG_ECM,EVf)

    wage   = nag.A_bar
    Vf_new = Array{Float64}( nag.nb, nag.nl, nag.nq, nag.nτ)
    bp_new = Array{Float64}( nag.nb, nag.nl, nag.nq, nag.nτ)
    cc_new = Array{Float64}( nag.nb, nag.nl, nag.nq, nag.nτ)

    # In what follows, we start at the end of the transitions, where B_{t+1}=B(τ1), and move backwards.
    for (i_τ, τ_v) in enumerate(nag.τgrid)
        for (i_q, q_v) in enumerate(nag.qgrid)

            for (i_l, l_v) in enumerate(nag.lgrid)
                τl = 1.0 - τ_v*l_v^(-nag.τ2)
                for (i_b, b_v) in enumerate(nag.bgrid)

                    rhs0 = -1e14
                    bpol = -1.0
                    cpol = 0.0

                    for (i_bp, bp_v) in enumerate(nag.bgrid)

                        cc    = wage*l_v*(1.0-nag.τY) + b_v  - τ_v - bp_v*q_v
                        if nag.prog_tax
                            cc    = wage*l_v*(1.0-τl) + b_v - bp_v*q_v
                        end
                        if cc <=0.0 ; rhs = -1e14
                        else
                            rhs   = ufun(nag,cc) + nag.β*EVf[i_bp, i_l, i_q, i_τ]
                        end

                        if rhs>rhs0
                            rhs0 = rhs
                            bpol = bp_v
                            cpol = cc
                        end

                    end

                    Vf_new[i_b,i_l,i_q,i_τ] = rhs0
                    bp_new[i_b,i_l,i_q,i_τ] = bpol
                    cc_new[i_b,i_l,i_q,i_τ] = cpol

                end
            end
        end
    end

    copy!(nag.Vf,Vf_new); copy!(nag.bprime_mat,bp_new); copy!(nag.c_mat,cc_new);

    maximum(bp_new)>0.0 || throw(error("bp too low at compute_bonds_in_grid"))

end


function iterate_VF_ongrid_nag!(nag::NAG_ECM)

    iter = 0
    dist = 10.0

    while iter<nag.maxit && dist>nag.tol_dist

        iter += 1

        Vf_mat_old = copy(nag.Vf)
        bp_old = copy(nag.bprime_mat)

        EVf = copy(nag.Vf)
        compute_EVf!(nag,EVf)
        compute_bonds_in_grid!(nag,EVf)

        dist   = maximum(abs,nag.Vf-Vf_mat_old)
        dist_b = maximum(abs,nag.bprime_mat-bp_old)

        if iter%50 == 0
            @printf("Iteration # %d for with distance in Vf of %.5f and in bonds of %.5f\n",iter,dist,dist_b)
        end
    end

    return Void
end
