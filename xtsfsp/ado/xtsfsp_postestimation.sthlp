{smcl}
{* *! version 1.0.0  10mar2015}{...}
{* revised: }{...}
{cmd:help xtsfsp postestimation}{right:also see: {help xtsfsp}}
{hline}

{title:Title}

{p2colset 5 32 38 2}{...}
{p2col :{hi:xtsfsp postestimation} {hline 2}}Postestimation tools for xtsfsp{p_end}
{p2colreset}{...}


{title:Description}

{pstd}
The following postestimation commands are available for {cmd:xtsfsp}.

{synoptset 13}{...}
{p2coldent :command}description{p_end}
{synoptline}
{synopt:{bf:{help estat}}}AIC, BIC, VCE, and estimation sample summary{p_end}
INCLUDE help post_estimates
INCLUDE help post_lincom
INCLUDE help post_lrtest
INCLUDE help post_nlcom
{synopt :{helpb sfsd postestimation##predict:predict}}predictions, residuals, influence statistics, and other diagnostic measures{p_end}
INCLUDE help post_predictnl
INCLUDE help post_test
INCLUDE help post_testnl
{synoptline}
{p2colreset}{...}


{marker predict}{...}
{title:Syntax for predict}

{p 8 16 2}{cmd:xtsfsp_p} {newvar} {ifin} [{cmd:,} {it:statistic}]

{synoptset 15 tabbed}{...}
{synopthdr :statistic}
{synoptline}
{syntab:Main}
{synopt :{opt xb}} prediction of the dependent variable; the default{p_end}
{synopt :{opt residuals}}estimates of residuals (depvar - xb) {p_end}
{synopt :{opt u}}estimates of (technical or cost) inefficiency via {it:E}(u|e) (Orea and Álvarez, 2019){p_end}
{synopt :{opt su}}estimates of spatial corrected inefficiency:[I-rho*W]^{-1}u{p_end}
{synopt :{opt te}}estimates of (technical or cost) efficiency via exp[-E(u|e)]{p_end}
{synopt :{opt ste}}estimates of spatial corrected efficiency{p_end}
{synopt :{opt uts}}estimates of commen inefficency u*_t{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
These statistics are only available for the estimation sample.




{title:Options for predict}

{dlgtab:Main}

{phang}
{opt xb}, the default, calculates the linear prediction.

{phang}
{opt residuals} calculates the composite residuals (depvar - xb).


{phang}
{opt u} produces estimates of (technical or cost) inefficiency via E(u|e) using the Orea and Álvarez (2019) estimator. 

{phang}
{opt su} produces spatial corrrected inefficiency for models with 
the spatial lag term. It is estimated by [I-rho*W]^{-1}E(u_t|e_t).
 
{phang}
{opt te} produces estimates of (technical or cost) efficiency via exp(-E(u|e)). 

{phang}
{opt ste} produces estimates of spatial corrected efficiency via  exp(-E(su|e)).

{title:Remarks}

{pstd} For the postestimation, the option {cmd:genwvars} should be specified in {cmd:xtsfsp} when estimating the models.{p_end}


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
{phang2}{bf:. {stata "xtsfsp y x, uhet(z) wu(w1,mata) wy(w1,mata) wv(w1,mata) wx(w1,mata) wxvars(x) genwvars"}}{p_end}

{phang2}{bf:. {stata "xtsfsp_p uhat, u"}}{p_end}

{phang2}{bf:. {stata "xtsfsp_p te, te"}}{p_end}
    


{marker author}{...}
{title:Author}

{pstd}
Kerui Du{break}
Xiamen University{break}
School of Management{break}
China{break}
{browse "kerrydu@xmu.edu.cn":kerrydu@xmu.edu.cn}{break}



{title:Also see}

{psee}
{space 2}Help:  {help xtsfsp} {help xtsfsp margins}.
{p_end}
