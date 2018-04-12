include("tauchen.jl")

using Interpolations, Optim, FastGaussQuadrature

type NAG_ECM

    # Model Parameters
    β::Float64  # Discount factor
    γ::Float64  # Risk aversion
    ρl::Float64 # Autocorrelation idiocyncratic shock
    σl::Float64 # std of idiocyncratic shock
    τY::Float64 # Tax on wages
    τ2::Float64 # Progressivity of taxation

    # Grid Parameters
    nl::Int  # Number of gridpoints for idiosyncratic shock
    nb::Int  # Number of gridpoints for bonds
    nq::Int  # Number of gridpoints for Price of Bonds
    nτ::Int  # Number of gridpoints for Lump-Sum Taxes
    nB::Int  # Number of gridpoints for cross-sectional mean
    nσB::Int # Number of gridpoints for cross-sectional s.d.
    nq0::Int # Index for inclusion
    num_nodes::Int # Number of points in histogram
    n_et::Int

    bgrid::Array{Float64, 1}  # Grid for bonds
    lgrid::Array{Float64, 1}  # Grid for idiosyncratic shock
    qgrid::Array{Float64, 1}  # Grid for bond prices in extended state space
    τgrid::Array{Float64, 1}  # Grid for lump-sum taxes
    Bgrid::Array{Float64, 1}  # Grid for cross-sectional mean
    σBgrid::Array{Float64, 1} # Grid for cross-sectional s.d.
    quad_nodes::Array{Float64}
    quad_weights::Array{Float64}

    # Finer grid for histogram
    b_grid_fine::Array{Float64}
    snodes::Array{Float64, 2}

    # Exogenous Objects
    A_bar::Float64 # mean TFP value
    L_bar::Float64 # mean Labor
    e_bar::Float64 # Coefficient on expenditure
    repay::Float64 # Debt repayment

    # Transition Matrices
    Πl::Array{Float64, 2} # Markov transition matrix for idiosyncratic shock
    Πl_ast::Array{Float64, 1} # Invariant distribution for idiosyncratic shock

    # Household's Value functions
    Vf::Array{Float64, 4}     # Value function
    Vb_mat::Array{Float64, 4} # Envelope Condition for Bonds

    # Household's Policies
    c_mat::Array{Float64, 4}
    bprime_mat::Array{Float64, 4}
    B_mat::Array{Float64, 2}

    # Options
    maxit::Int64
    tol_dist::Float64

end

function nag_economy(;  β   = 0.97,
                        γ   = 2.0,
                        ρl  = 0.9,
                        σl  = 0.05,
                        τY  = 0.005,
                        nl  = 7,
                        nb  = 151,
                        nq  = 20,
                        nσB = 5,
                        nq0 = nq,
                        nτ  = 15, #9,
                        nB  = nτ,
                        n_et = 10,
                        A_bar = 1.0,
                        L_bar = 1.0,
                        e_bar = 0.25+τY, #0.01+τY,
                        repay = 1.0,
                        maxit = 500,#2000,
                        tol_dist = 1e-4,#5e-3
                        τ2 = 0.10)


    # Calibration
    # ρl    = 0.6   # Aiyagari (1994)
    # σl    = 0.2   # Aiyagari (1994)

    do_MCG = 0

    if do_MCG == 1
        τY    = 0.376 # Aiyagari - McGrattan (1998) for b = 2/3 and r = 4.5%
        e_bar = 0.217 # Aiyagari - McGrattan (1998) for b = 2/3 and r = 4.5%
    end


    """ Define grids """
    bgrid     = collect(linspace(1e-3,9.0,nb))
    qgrid     = collect(linspace(β+0.001,1.0/1.001,nq))
    Bgrid     = collect(linspace(0.01,3.0,nB))
    σBgrid    = collect(linspace(0.005, 0.5, nσB))

    quad_nodes, quad_weights = gausshermite(n_et)


    """ Define process for endowment of labor """
    lgrid, Πl = tauchen_fun(nl, ρl, σl)
    lgrid     = collect(exp.(lgrid))

    # # Based on Davila et al (2012)
    # lgrid   = [.65;1.0;1.7]
    # Πl      = [[0.99 1.0-0.99 0.0]; [0.01 0.985 1-0.985-0.01]; [0.0 1-0.85 0.85]]

    # _, Πl_inv   = eigs(Πl', nev=1)
    # Πl_inv      = real(squeeze(Πl_inv, 2))
    # Πl_inv      = Πl_inv ./ sum(Πl_inv)

    # lgrid[2]    = (1.0 - Πl_inv[1]*lgrid[1] - Πl_inv[3]*lgrid[3])/Πl_inv[2]
    # nl          = length(lgrid)

    _, Πl_ast   = eigs(Πl', nev=1)
    Πl_ast      = real(squeeze(Πl_ast, 2))
    Πl_ast      = Πl_ast ./ sum(Πl_ast)

    τmin1 = e_bar-τY+1e-8
    τmax1 = A_bar*lgrid[1]*(1.0-τY)*0.8
    τgrid = collect(linspace(τmin1,min(τmax1,0.3),nτ))
    # τgrid = collect(linspace(τmin1,τmin1+0.01,nτ))

    if do_MCG == 1
        τgrid = collect(linspace(-0.145,-0.09,nτ))
    end

    # Initial guess for Envelope
    Vb_mat = zeros(nb,nl,nq,nτ)

    for (i_τ,τ_v) in enumerate(τgrid)
        for (i_q,q_v) in enumerate(qgrid)
            for (i_l,l_v) in enumerate(lgrid)
                τl = 1.0 - τ_v*l_v^(-τ2)
                for (i_b,b_v) in enumerate(bgrid)
                    disp_income = A_bar*l_v*(1.0-τY) - τ_v + b_v*repay
                    disp_income = A_bar*l_v*(1.0-τl) + b_v*repay
                    disp_income>0.0 || throw(error("Negative disp_income"))
                    Vb_mat[i_b,i_l,i_q,i_τ] = (disp_income*0.9)^(-γ)
                end
            end
        end
    end

    bprime_mat,c_mat,Vf = zeros(Vb_mat),zeros(Vb_mat),zeros(Vb_mat)

    b_grid_fine = collect(linspace(bgrid[1],bgrid[end],300))
    snodes      = [repmat(b_grid_fine,length(lgrid),1) kron(lgrid,ones(length(b_grid_fine)))]

    num_nodes = size(snodes,1)

    B_mat = Array{Float64}(nq,nτ)

    for i_τ=1:nτ
        for i_q=1:nq
            B_mat[i_q,i_τ] = (e_bar - τY - τgrid[i_τ])/(qgrid[i_q] - 1.0)
            tax_coll = 1.0 - τ_v * dot(lgrid.^(1.0-τ2), Πl_ast)
            B_mat[i_q,i_τ] = (e_bar - tax_coll)/(qgrid[i_q] - 1.0)
        end
    end

    return NAG_ECM(β,γ,ρl,σl,τY,nl,nb,nq,nτ,nB,nσB,nq0,num_nodes,n_et,bgrid,lgrid,qgrid,τgrid,Bgrid,σBgrid,quad_nodes,quad_weights,b_grid_fine,snodes,
                   A_bar,L_bar,e_bar,repay,Πl,Πl_ast,Vf,Vb_mat,c_mat,bprime_mat,B_mat,maxit,tol_dist)

end


e_fun(nag::NAG_ECM,A_v)       = nag.e_bar*A_v
ufun(nag::NAG_ECM,cons)       = cons.^(1.0-nag.γ)./(1.0-nag.γ)
uprime(nag::NAG_ECM,cons)     = cons.^(-nag.γ)
uprime_inv(nag::NAG_ECM,Vval) = Vval.^(-1.0/nag.γ)
u_inv(nag::NAG_ECM, V)        = (V*(1.0-nag.γ)).^(1.0/(1.0-nag.γ))
cons_eq(nag::NAG_ECM, V)      = u_inv(nag, (1.0-nag.β)*V)

# """ Solve the Household's Problem """

include("households_problem_ssnag_vq.jl")

nag        = nag_economy()
ctilde_mat = zeros(size(nag.Vb_mat))
bprime_mat = copy(ctilde_mat)

@time iterate_envelope_nag!(nag, ctilde_mat)
@time compute_consumption_nag!(nag)
@time compute_Vf_nag!(nag)


numerator, denominator = 0.0, 0.0

for ii=1:nag.nl, jj=1:nag.nl
    numerator   += abs(nag.lgrid[ii]-nag.lgrid[jj])
    denominator += nag.lgrid[jj]
end

gini_index_labor = numerator / (2*denominator)



# """ Solve for the SS equilibrium """

include("solve_ssnag_eqm_vq.jl")

# τvec_comp_an = [nag.τgrid[1];nag.τgrid[end]]
τvec_comp_an = copy(nag.τgrid)

function compute_all_SS(nag::NAG_ECM, τvec_comp_an; do_simple_in=true)

    nτ       = size(τvec_comp_an,1)
    λast_out = Array{Float64}(size(nag.snodes,1),nτ)
    dens_b   = Array{Float64}(length(nag.b_grid_fine),nτ)
    B_demand, B_supply, q_out_SS, welf_out = Array{Float64}(nτ), Array{Float64}(nτ), Array{Float64}(nτ), Array{Float64}(nτ)

    qgrid_cell_SS, bp_cell_SS, Vf_cell_SS, Vb_cell_SS = [], [], [], []

    for i_τ=1:nτ

        out_SS_objects  = compute_ss_root(nag, τvec_comp_an[i_τ], do_simple=do_simple_in)

        λast_out[:,i_τ] = out_SS_objects[1]
        B_demand[i_τ]   = out_SS_objects[2]
        B_supply[i_τ]   = out_SS_objects[3]
        q_out_SS[i_τ]   = out_SS_objects[4]
        qgrid_input     = out_SS_objects[5]
        Vf_input        = out_SS_objects[6]
        bpmat_out_SS    = out_SS_objects[7]
        Vbp_out_SS      = out_SS_objects[8]

        welf_out[i_τ] = compute_welfare_nag(nag, λast_out[:,i_τ], q_out_SS[i_τ], qgrid_input, Vf_input)
        dens_b[:,i_τ] = sum(reshape(λast_out[:,i_τ],length(nag.b_grid_fine),nag.nl),2)

        push!(qgrid_cell_SS, qgrid_input)
        push!(bp_cell_SS,    bpmat_out_SS)
        push!(Vf_cell_SS,    Vf_input)
        push!(Vb_cell_SS,    Vbp_out_SS)


    end

    return λast_out, dens_b, B_demand, B_supply, q_out_SS, welf_out, qgrid_cell_SS, bp_cell_SS, Vf_cell_SS, Vb_cell_SS
end

@time λast_out, dens_b, B_demand, B_supply, q_out_SS, welf_out, qgrid_cell_SS, bp_cell_SS, Vf_cell_SS, Vb_cell_SS = compute_all_SS(nag, τvec_comp_an; do_simple_in=true);

nag.Bgrid = copy(B_demand)

# # ## Compute Gini Indices on Wealth
# # income_dens = zeros(nag.num_nodes,3)

# # for i_s = 1:nag.num_nodes
# #     income_dens[i_s, 1] = nag.snodes[i_s,1] + nag.snodes[i_s,2]
# #     income_dens[i_s, 2] = λast_out[i_s, 1]
# #     income_dens[i_s, 3] = λast_out[i_s, 2]
# # end

# # income_dens_sorted = sortrows(income_dens)

# # numerator1, numerator2 = 0.0, 0.0
# # for i_s = 1:nag.num_nodes

# #     Sim1_1, S_i_1 = 0.0, 0.0
# #     Sim1_2, S_i_2 = 0.0, 0.0

# #     if i_s > 1
# #         for jj=1:i_s-1
# #             Sim1_1 += income_dens_sorted[jj, 2]*income_dens_sorted[jj, 1]
# #             Sim1_2 += income_dens_sorted[jj, 3]*income_dens_sorted[jj, 1]
# #         end
# #     end
# #     for jj=1:i_s
# #         S_i_1  += income_dens_sorted[jj, 2]*income_dens_sorted[jj, 1]
# #         S_i_2  += income_dens_sorted[jj, 3]*income_dens_sorted[jj, 1]
# #     end

# #     numerator1 += income_dens_sorted[i_s, 2]*(Sim1_1 + S_i_1)
# #     numerator2 += income_dens_sorted[i_s, 3]*(Sim1_2 + S_i_2)
# # end

# # denominator1, denominator2 = 0.0, 0.0

# # for jj=1:nag.num_nodes
# #     denominator1  += income_dens_sorted[jj, 2]*income_dens_sorted[jj, 1]
# #     denominator2  += income_dens_sorted[jj, 3]*income_dens_sorted[jj, 1]
# # end

# # gini_wealth     = zeros(2)

# # gini_wealth[1]  = 1.0 - numerator1/denominator1
# # gini_wealth[2]  = 1.0 - numerator2/denominator2

# # @printf("Gini indices on wealth are %.2f and %.2f\n", gini_wealth[1], gini_wealth[2])


# # b_poster = [0.0;1.0;3.0]


# # include("do_figures_ssnag_vq.jl")
# # include("compute_haircut.jl")
# # include("do_figures_ssnag_trandyn.jl")

