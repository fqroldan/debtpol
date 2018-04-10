# include("ssnag_main_vq.jl")

using PyPlot

fig, axis = subplots(1, 3, figsize = (10,5))
ax = axis[1]; ax[:plot](100*B_supply, welf_out);ax[:set_xlabel]("B/Y")
ax[:set_title]("Welfare")
ax = axis[2]; ax[:plot](100*B_supply, (1./q_out_SS-1.0)*100); ax[:set_xlabel]("B/Y")
ax[:set_ylabel](L"\%"); ax[:set_title]("Interest Rates")
ax = axis[3]; ax[:plot](100*B_supply, τvec_comp_an*100); ax[:set_xlabel]("B/Y")
ax[:set_ylabel]("% of GDP"); ax[:set_title]("Taxes")
tight_layout()
# savefig(pwd() * "/../../../output/steady_state/welfare_new.png")
# savefig(pwd() * "/../../../output/steady_state/welfare_new.pdf")


fig, ax = subplots()
ax[:scatter](100*B_supply, B_demand-B_supply, label = L"B^d-B^s"); ax[:set_xlabel](L"B")
ax[:set_title]("Excess demand"); ax[:legend](loc=1)
tight_layout()
# savefig(pwd() * "/../../../output/steady_state/excess_dem_at_welfare_new.png")
