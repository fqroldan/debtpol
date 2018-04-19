# include("ssnag_main_vq.jl")

usar_plotly = true

if usar_plotly
	using PlotlyJS, Rsvg
	p1 = plot(scatter(; x=100*B_supply, y=welf_out, showlegend=false), 			  Layout(; title="Welfare", 	   xaxis=attr(title="B/Y")));
	p2 = plot(scatter(; x=100*B_supply, y=(1./q_out_SS-1.0)*100, showlegend=false), Layout(; title="Interest Rates", xaxis=attr(title="B/Y"), yaxis=attr(title="%")));
	p3 = plot(scatter(; x=100*B_supply, y=τvec_comp_an*100, showlegend=false), 	  Layout(; title="Taxes", 		   xaxis=attr(title="B/Y"), yaxis=attr(title="% of GDP")));

	p = [p1 p2 p3]
	p.plot.layout["width"] = 800
	p.plot.layout["height"] = 400
	p.plot.layout["font_family"] = "Fira Sans Light"
	savefig(p, pwd() * "/../../../output/steady_state/welfare_new_prog.pdf")
	p
else
	using PyPlot, LaTeXStrings

	fig, axis = subplots(1, 3, figsize = (10,5))
	ax = axis[1]; ax[:plot](100*B_supply, welf_out);ax[:set_xlabel]("B/Y")
	ax[:set_title]("Welfare")
	ax = axis[2]; ax[:plot](100*B_supply, (1./q_out_SS-1.0)*100); ax[:set_xlabel]("B/Y")
	ax[:set_ylabel](L"\%"); ax[:set_title]("Interest Rates")
	ax = axis[3]; ax[:plot](100*B_supply, τvec_comp_an*100); ax[:set_xlabel]("B/Y")
	ax[:set_ylabel]("% of GDP"); ax[:set_title]("Taxes")
	tight_layout()
	savefig(pwd() * "/../../../output/steady_state/welfare_new_prog.png")


	fig, ax = subplots()
	ax[:scatter](100*B_supply, B_demand-B_supply, label = L"B^d-B^s"); ax[:set_xlabel](L"B")
	ax[:set_title]("Excess demand"); ax[:legend](loc=1)
	tight_layout()
	savefig(pwd() * "/../../../output/steady_state/excess_dem_at_welfare_new_prog.png")
end

include("do_figures_ssnag_vq.jl")