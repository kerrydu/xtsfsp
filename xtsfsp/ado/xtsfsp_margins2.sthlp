{smcl}
{cmd:help xtsfsp_margins}
{hline}

{title:Title}

{p2colset 5 13 15 2}{...}
{p2col :{hi:xtsfsp_margins} {hline 2}}compute marginal effects 
of Spatial panel stochastic frontier models in the style of {help xtsfsp##OA2019:{bind:Orea and Álvarez (2019)}} and {help xtsfsp##Galli2022:{bind:Galli (2022)}} {p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{pstd}
Syntax

{p 8 17 2}
{cmd:xtsfsp_margins} {varname} [{cmd:,} {it:options}]


{synoptset 31 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt mc(#)}}set # of replications for Monte Carlo simulation{p_end}
{synopt :{opt seed(#)}}set seed for generating random number{p_end}
{synopt :{nodots}}suppress iteration dots{p_end}



{marker description}{...}
{title:Description}

{pstd}
{opt xtsfsp_margins} computes the total, direct and indirect marginal effects for spatial panel stochastic production or cost frontier models
proposed by
{help xtsfsp##OA2019:{bind:Orea and Álvarez (2019)}} and {help xtsfsp##Galli2022:{bind:Galli (2022)}}. See 
{help xtsfsp##OA2019:{bind:Orea and Álvarez (2019)}} and {help xtsfsp##Galli2022:{bind:Galli (2022)}} for a detailed
explanation of their methodology and empirical analyses.{p_end}
{p2colreset}{...}


{marker options}{...}
{title:Options}

{phang}
{cmd:mc(#)} specifes # of replications for Monte Carlo simulation
for estimating the standard errors. If not specified, the delta method is used.

{phang}
{cmd:seed(#)} specifies the seed for Monte Carlo simulation. The default is 123.

{phang}
{cmd: nodots} suppress iteration dots.

{marker examples}{...}
{title:Examples}

    {title:yxuv-SAR model with time-invariant spatial weight matrix}

{pstd}
Setup{p_end}
{phang2}{bf:. {stata "mata mata matuse xtsfsp_w1,replace"}}{p_end}
{phang2}{bf:. {stata "use xtsfsp_ex1.dta"}}{p_end}
{phang2}{bf:. {stata "xtset id t"}}{p_end}

{pstd}
Stochastic production model with four different sources of spatial cross-sectional dependence {p_end}
{phang2}{bf:. {stata "xtsfsp y x, uhet(z) wu(w1,mata) wy(w1,mata) wv(w1,mata) wx(w1,mata) wxvars(x) te(te)"}}{p_end}

{phang2}{bf:. {stata "xtsfsp_margins z"}}{p_end}

{phang2}{bf:. {stata "xtsfsp_margins z, mc(200)"}}{p_end}
    
{marker results}{...} 
{title:Stored results}

{pstd}
{cmd:xtsfsp_margins} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(totalxmargins)}}table of total marginal effect for variabls in frontier{p_end}
{synopt:{cmd:e(directxmargins)}}table of direct marginal effect for variabls in frontier{p_end}
{synopt:{cmd:e(indirectxmargins)}}table of indirect marginal effect for variabls in frontier{p_end}
{synopt:{cmd:e(totalzmargins)}}table of total marginal effect for variabls in scaling function{p_end}
{synopt:{cmd:e(directzmargins)}}table of direct marginal effect for variabls in scaling function{p_end}
{synopt:{cmd:e(indirectzmargins)}}table of indirect marginal effect for variabls in scaling function{p_end}

{marker acknowledgments}{...}
{title:Acknowledgments}

{pstd}
Kerui Du acknowledges financial support from the National Natural Science Foundation of China(Grant no. 72074184).



{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
{cmd:xtsfsp} is not an official Stata command.  It is a third-party command
programmed as a free contribution
to the research society.  By choosing to download, install, and use the
{cmd:xtsfsp} package, users assume all the liability for any
{cmd:xtsfsp}-package-related risk.  If you encounter any problems with the
{cmd:xtsfsp} package, or if you have comments, suggestions, or questions, please
send an email to Kerui Du at 
{browse "mailto:kerrydu@xmu.edu.cn":kerrydu@xmu.edu.cn}.


{marker citation}{...}
{title:References}

{marker OA2019}{...}
{phang}
Orea, Luis, and Inmaculada C. Álvarez. “A New Stochastic Frontier Model with Cross-Sectional Effects in Both Noise and Inefficiency Terms.” 
Journal of Econometrics 213, 
no. 2 (December 1, 2019): 556–77. https://doi.org/10.1016/j.jeconom.2019.07.004.

{marker Galli2022}{...}
{phang}
Galli, Federica. “A Spatial Stochastic Frontier Model Including Both Frontier and Error-Based Spatial 
Cross-Sectional Dependence.” 
Spatial Economic Analysis 0, no. 0 (July 26, 2022): 1–20. 
https://doi.org/10.1080/17421772.2022.2097729.

{phang}
Belotti, F., Daidone, S., Ilardi, G.,  Atella, V. "Stochastic frontier analysis using Stata", Stata Journal, 2013, 13(4):719-758.

{phang}
 Karakaplan, M. U. "Estimating endogenous stochastic frontier models in Stata", Stata Journal, 2017,17: 39-55. 

{phang}
Ditzen, Grieser, Zekhnini. "nwxtregress - network regression in Stata", 2023. https://github.com/JanDitzen/nwxtregress



{marker author}{...}
{title:Author}

{pstd}
Kerui Du{break}
Xiamen University{break}
School of Management{break}
China{break}
{browse "kerrydu@xmu.edu.cn":kerrydu@xmu.edu.cn}{break}


{pstd}
Luis Orea{break}
University of Oviedo{break}
Department of Economics{break}
Spain{break}
{browse "lorea@uniovi.es":lorea@uniovi.es}{break}


{pstd}
Inmaculada C. Álvarez{break}
Universidad Autónoma de Madrid{break}
Department of Economics{break}
Spain{break}
{browse "inmaculada.alvarez@uam.es":inmaculada.alvarez@uam.es}{break}


{marker see}{...}
{title:Also see}

{p 7 14 2}{manhelp frontier R}, 
{manhelp xtfrontier XT}{p_end} 
{p 7 14 2}{helpb sfpanel}, {helpb sfkk},  {helpb nwxtregress} (if installed)
