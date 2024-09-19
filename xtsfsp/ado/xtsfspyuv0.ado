* version 1.0, 2023-10-10
* subprogram for xtsfsp 
* Case8: y=rhoWy+x+v-u; v=[I-gammaW]^{-1}v*; u=[I-tauW]^{-1}u*
capture program drop xtsfspyuv0
program define xtsfspyuv0, eclass sortpreserve
version 16

// 	if replay() {
// 		if (`"`e(cmd)'"' != "xtsfspyuv0") error 301
// 		Replay `0'
// 	}
// 	else	Estimate `0'
// end

// program Estimate, eclass sortpreserve

syntax varlist,  Uhet(string) [ INItial(name) NOCONstant NORMalize(string) wu(string) ///
                              wv(string) te(name) GENWVARS mldisplay(string) ///
                              DELmissing MLPLOT NOGraph MLMODELopt(string)  ///
							  level(real 95) COST wxvars(varlist) genwxvars ///
							  MLSEarch(string) MLMAXopt(string) DELVE ///
							  CONSTraints(string) wy(string) wx(string) ///
							  LNDETFULL lndetmc(numlist >0 min=2 max=2) NOLOG Vhet(string)] 

local varlist: list uniq varlist
local uhet: list uniq uhet
local vhet: list uniq vhet
local wxvars: list uniq wxvars
if ("`nolog'"!="") local nolog qui
//marksample touse 
local cmdline xtsfsp `0'
if "`genwxvars'" !="" local genwvars genwvars
preserve
marksample touse 

local diopts level(`level') `mldisplay'
mlopts std, `mlmodelopt' `mlmaxopt' `constraints'
local cns constraints(`constraints')
//gettoken yvar xvars: varlist 
//_fv_check_depvar `yvar'

if ("`wxvars'"=="" & "`wx'"!=""){
	//di  `"varlist is not specified in wxvars(), wx(`wx') is neglected"'
	di as error `"wx(`wx') must be combined with wxvars()"'
	exit 198
}  

if ("`wxvars'"!="" & "`wx'"!="" & "`genwvars'"!=""){
	foreach xv in `xvars'{
		confirm new var W_`xv'
	}
}
gettoken yvar xvars: varlist
if ( "`genwvars'"!=""){
	confirm new var W_`yvar'
}

if ("`initial'"!="" & "`delve'"!=""){
	di "Warning: initial(`initial') overrides delve"
}
if ("`initial'"!="" & "`mlsearch'"!=""){
	di "Warning: initial(`initial') overrides mlsearch(`mlsearch')"
}
if ("`delve'"!="" & "`mlsearch'"!=""){
	di "Warning: delve overrides mlsearch(`mlsearch')"
}



parsespmat0 `wu' 
parsespmat1 `wu' `r(ldot)' aname(w_ina)
local nwu = r(nw)

// check data with spmatrix
	_xt, trequired 
	local id=r(ivar)
	local time=r(tvar)
	gettoken uhet0 : uhet, p(,)
	gettoken vhet0 : vhet, p(,)
	markout `touse' `uhet0' `vhet0'
    qui keep `varlist' `wxvars' `id' `time' `uhet0' `touse' `vhet0'
    tempvar order0
    qui gen int `order0' =_n
// sort data	
	qui issorted `time' `id'	
	//tempvar time2

	qui distinct2 `time'
	local T = r(ndistinct)	
//mata mata describe
//	mata: marksuse = st_data(.,"`touse'")
	qui count if `touse'==0
	local nummissing = r(N)
	if(`nummissing'>0){
		mata: marksuse = st_data(.,"`touse'")
	}    
	mata: _order_0   = st_data(.,"`order0'","`touse'") // record the row# in the original data	

	if (`nwu'!=1 & `nwu'!=`T') {
		di as error "Spatial weight matrixs in wu() are specified as time-varying, but # of spmatrix != # of periods"
		exit 198
	}    
	tempvar time2
	qui egen `time2' = group(`time')
	global paneltvar `time2'
	//mata mata describe
	checkspmat w_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')

    scalar rmin = max(-0.9999,r(min_w_ina))
	scalar rmax = min(0.9999,r(max_w_ina))
	global rmin = rmin
	global rmax = rmax

	if ("`wx'"=="`wy'"){
		local wxwx w_ina
	}
	else{
		local wxwx wx_ina
	}
    if ("`wx'"!="" & `"`wxvars'"'!=""  & "`wxwx'"=="wx_ina"){
        parsespmat0 `wx' 
        parsespmat1 `wx' `r(ldot)' aname(wx_ina)
        local nwx = r(nw)
        if (`nwx'!=1 & `nwx'!=`T') {
            di as error "Spatial weight matrixs in wx() are specified as time-varying, but # of spmatrix != # of periods"
            exit 198
        }  
        qui checkspmat wx_ina, time(`time2') touse(`touse')  `delmissing' normalize(`normalize')
    }



    qui keep if `touse'
    mata: _pan_tvar =st_data( .,"`time2'")	
	
	
	 	if ("`lndetfull'"!=""){
			local bp bp 
			mata: _rho_lndet_ = panlndetfull(w_ina,$rmin,$rmax,`T')
		}
		if ("`lndetmc'"!=""){
			local bp bp 
			tokenize `lndetmc'
			mata: _rho_lndet_ = panlndetmc(`1',`2',w_ina,$rmin,$rmax,`T')
		}	

**************
	qui genwxvars `yvar', aname(w_ina) tvar(`time2') // Wednesday, June 26, 2024 at 21:35:25
	local wyvar  `r(wxnames)'
	mata: _order_wyvar = st_data(.,"`wyvar'","`touse'")	
   * generating Wx
	if(`"`wxvars'"'!=""){
      qui genwxvars `wxvars', aname(`wxwx') tvar(`time2')
      local wxvars2  `r(wxnames)'
      mata: _order_wx = st_data(.,"`wxvars2'","`touse'")
	  //cap mata mata drop wx_ina
	}	
************

	if("`initial'"=="" & "`delve'"!="") { 
		//qui genwxvars `yvar', aname(w_ina) tvar(`time2')
	    //local wyvar `r(wxnames)'
		qui corr `yvar' `wyvar'
		local rhoy = r(rho)
		qui frontier `yvar' `xvars' `wxvars2',`noconstant' uhet(`uhet') iterate(50) `cns' vhet(`vhet')
	    mat b0 =e(b)
		qui predict double xbfron
		qui gen double ehfron = `yvar'-xbfron
		qui genwxvars ehfron, aname(w_ina) tvar(`time2')
	    local wehfron `r(wxnames)'
		qui corr ehfron `wehfron'
		local rhouv = r(rho)
		//local r0 = runiform(rmin,rmax)
		//local r0 = 0.3
		//local r1 = (`r0'-rmin)/(rmax-rmin)
		//mat b0=b0, ln(`r1'/(1-`r1')),ln(`r1'/(1-`r1')),ln(`r1'/(1-`r1'))
		mat b0=b0, ln(`rhoy'/(1-`rhoy')),ln(`rhouv'/(1-`rhouv')),ln(`rhouv'/(1-`rhouv'))
	}

	if ("`cost'"!=""){
		mata: _cost = -1
	} 
	else{
		mata: _cost = 1
	}	
    if `"`vhet'"'=="" {
		local vterm /lnsigv2
	}
	else{
		local vterm (lnsigv2: `vhet') 
	}
	local modeltype = cond("`wxvars'"=="","yuv-SAR","yxuv-SAR")
	local title Spatial frontier model(`modeltype')

// eq1: frontier; eq2: ln(sigmav2); eq3: uhet; eq4:ln(sigmau2); eq5: Wy; eq6: Wv; eq7:Wu	
	ml model d0 xtsfspyuv0`bp'() (frontier:`yvar' = `xvars' `wxvars2',`noconstant') ///
	              `vterm'   (uhet: `uhet') ///
	            (Wy:) (Wv:) (Wu:), nopreserve `cns' `mlmodelopt'  title(`title')
	
	
	if("`initial'"=="" & "`delve'"!="") { 
		ml init b0,copy
		//ml init Wu:_cons=0.2
		//ml init Wv:_cons=0.2
	}
	if ("`initial'"=="" & "`delve'"==""){
		`nolog' display "searching initial values..."
		`nolog' ml search, `mlsearch'
	} 
	if ("`initial'"!="") ml init `initial', copy
	if ("`mlplot'"!=""){
		`nolog' display "mlplot initial values..."
		if "`nograph'"!="" set graphics off
		`nolog'  ml plot uhet:_cons
		`nolog' cap ml plot lnsigv2:_cons
		`nolog' cap ml plot /lnsigv2
		`nolog' ml plot Wy:_cons
		`nolog' ml plot Wu:_cons
		`nolog' ml plot Wv:_cons
		if "`nograph'"!="" set graphics on
	}
  	local mlmaxopt `mlmaxopt' noout difficult
	local mlmaxopt: list uniq mlmaxopt 
   `nolog' ml max, `mlmaxopt'  

   ereturn local cmd xtsfsp
   ereturn local cmdbase ml
   ereturn local cmdline `cmdline'
   ereturn local wy w_ina 
   ereturn local wv w_ina 
   ereturn local wu w_ina 
   ereturn local wx `wxwx'
   ereturn local wxvars `wxvars'
   ereturn local depvar `yvar' 
   ereturn local xeq `xvars' `wxvars2',`noconstant'
   ereturn local veq `vhet'
   ereturn local ueq `uhet'
   ereturn local hasgenwvars `genwvars'
   ereturn local  predict  "xtsfsp_p"
   ereturn local  margins  "xtsfsp_margins"   
   ereturn scalar T = `T'
   ereturn scalar rumin = $rmin
   ereturn scalar rumax = $rmax
   ereturn scalar rymin = $rmin
   ereturn scalar rymax = $rmax   
   ereturn scalar rvmin = $rmin
   ereturn scalar rvmax = $rmax  
   ereturn local ivar `id'
   ereturn local tvar `time'  
   local tau = _b[Wu:_cons]
   local tau = $rmin/(1+exp(`tau'))+$rmax*exp(`tau')/(1+exp(`tau'))
   local gamma = _b[Wv:_cons]
   local gamma = $rmin/(1+exp(`gamma'))+$rmax*exp(`gamma')/(1+exp(`gamma'))  
   local rho = _b[Wy:_cons]
   local rho = $rmin/(1+exp(`rho'))+$rmax*exp(`rho')/(1+exp(`rho'))
  
   ereturn scalar rho = `rho'
   ereturn scalar gamma = `gamma'
   ereturn scalar tau = `tau'


   Replay , `diopts'
   if `"`wxvars'"'!="" di "      W_(`wxvars') represent Spatial Durbin terms W(`wxvars')"

/*
   if(`"`te'"'!=""){
		tempname bml
		mat `bml' = e(b)
		mata: _b_ml = st_matrix("`bml'")	
	    local nx: word count `xvars' `wxvars2'
		//local nz: word count `uhet'
		if("`noconstant'"=="") local noconstant constant
		mata:_te_order=xtsfspyuv0_te(_b_ml,`nx',"`yvar'","`xvars' `wxvars2'","`uhet'","`vhet'","`noconstant'")
		mata: rho = $rmin/(1+exp(_b[Wy:_cons]))+$rmax*exp(_b[Wy:_cons])/(1+exp(_b[Wy:_cons]))	
		mata:spte(w_ina,rho,`T',_te_order,htildeut=.)
   }
*/

// if(`"`te'"'!=""){
// 	local wu w_ina 
// 	local wy w_ina
// 	local wv w_ina
// 	tempname bml
// 	mat `bml' = e(b)
// 	mata: _b_ml = st_matrix("`bml'")
// 	tempname yyy xxx vvv zzz
// 	data2mata `yyy' = `yvar',noconstant
// 	data2mata `xxx' = `xvars' `wxvars2',`noconstant'
// 	local nx = r(ncol)
// 	data2mata `vvv' = `vhet'
// 	local nv = r(ncol)
// 	data2mata `zzz' = `uhet'
// 	local nu = r(ncol)
// 	mata: bx = _b_ml[1..`nx'],`rho' // - rho*wy
// 	mata: `xxx' = `xxx',_order_wyvar // - rho*wy
// 	mata: bv = _b_ml[(`nx'+1)..(`nx'+`nv')]
// 	mata: bu = _b_ml[(`nx'+`nv'+1)..(`nx'+`nv'+`nu')]
// 	mata: _te_order = spte(`xxx',bx',`vvv',bv',`zzz',bu',`rho',`gamma',`tau',`yyy',`wy',`wv',`wu')
// }


  restore
  
  	if(`"`wxvars'"'!=""&"`genwvars'"!=""){
      foreach v in `wxvars'{
        qui gen double W_`v' = .
        label var W_`v' `"W*`v'"'
      }
	  mata: getdatafmata(_order_wx,_order_0,"`wxvars2'")
      cap mata mata drop  _order_wx
	  ereturn local wx `wxwx'

	}
	if("`genwvars'"!=""){
		qui gen double W_`yvar' = .
		label var W_`yvar' "W*`yvar'"
		mata: getdatafmata(_order_wyvar,_order_0,"W_`yvar'")
		//cap mata mata drop  _order_wyvar
	  }	
//    if(`"`te'"'!=""){
// 		qui gen double `te' = .
// 		label var `te' "technical efficiency"
// 		mata: getdatafmata(_te_order,_order_0,"`te'")
// 		cap mata mata drop  _te_order		
//    }  	
   	if(`nummissing'>0){
		di "      Missing values found"
		di "      The regression sample recorded by variable __e_sample__"
		cap drop __e_sample__
		qui cap gen byte __e_sample__ = 0
		label var __e_sample__ "e(sample)"
		mata: getdatafmata(J(length(_order_0),1,1),_order_0,"__e_sample__")
		//cap mata mata drop  _touse		
	}	 
   cap mata mata drop _order_0

end


///////////////////////subprograms//////////////////
cap program drop Replay
program Replay
	syntax [, Level(cilevel) *]
	ml display , level(`level')	`options'	///
		diparm(Wy, label(rho) prob function($rmin/(1+exp(@))+$rmax*exp(@)/(1+exp(@))) /*
       */ d(exp(@)*(($rmax-$rmin)/(1+exp(@))^2))) ///        
		diparm(Wv, label(gamma) prob function($rmin/(1+exp(@))+$rmax*exp(@)/(1+exp(@))) /*
       */ d(exp(@)*(($rmax-$rmin)/(1+exp(@))^2))) ///
	    diparm(Wu, label(tau) prob function($rmin/(1+exp(@))+$rmax*exp(@)/(1+exp(@))) /*
       */ d(exp(@)*(($rmax-$rmin)/(1+exp(@))^2)))
	   global diparmopt diparm(Wy, label(rho) prob function($rmin/(1+exp(@))+$rmax*exp(@)/(1+exp(@)))  d(exp(@)*(($rmax-$rmin)/(1+exp(@))^2))) diparm(Wv, label(gamma) prob function($rmin/(1+exp(@))+$rmax*exp(@)/(1+exp(@)))  d(exp(@)*(($rmax-$rmin)/(1+exp(@))^2))) diparm(Wu, label(tau) prob function($rmin/(1+exp(@))+$rmax*exp(@)/(1+exp(@)))  d(exp(@)*(($rmax-$rmin)/(1+exp(@))^2)))
	   global end1 "Note: Wy:_cons, Wv:_cons and Wu:_cons are the transformed parameters;"
	   global end2 "      rho, gamma and tau are their origin metrics, respectively"
	   di "${end1}"
	   di "${end2}"
end


//include spxtsfa_aug.ado

//////utility comands and function for spxtsfa////

cap program drop genwxvars
program define genwxvars,rclass

version 16

syntax varlist, aname(string) [tvar(varname)]

if `"`tvar'"'==""{
	tempvar tvar 
	qui gen  byte `tavr'=1
}


mata: _genwxvars("`varlist'",`aname',"`tvar'")
return local wxnames `wxnames'

end


capture program drop checkspmat
program define checkspmat, rclass

syntax namelist(name=wnames), time(varname) touse(varname) [DELMissing NORMalize(string)]

//preserve

qui count if `touse'==0
local n0 = r(N)

if `n0'>0 & "`delmissing'"==""{
	di as red "missing values found. use delmissing to remove the units from the spmatrix"
	error 198
}


if `n0'>0 & "`delmissing'"!=""{
	di  "missing values found. The corresponding units are deleted from the spmatrix" _n
}

if "`normalize'"==""{
	local ntype=0
}
else if "`normalize'"=="row"{
	local ntype=1
}
else if "`normalize'"=="col"{
	local ntype=2
}
else if "`normalize'"=="spectral"{
	local ntype=3
}
else if "`normalize'"=="minmax"{
	local ntype=4
}
else{
	di as error "errors in normalize(), one of {row,col,spectral,minmax} should be specified. "
	error 198
}


//mata mata describe
foreach w in `wnames'{
    mata: _checkspmat("`time' `touse'",`w',`ntype')
    return scalar rmin_`w' = r(rmin)
    return scalar rmax_`w' = r(rmax)
}


end


///////////////////////////////////

capture program drop parsespmat0
program define parsespmat0,rclass
syntax namelist(name=wnames),[MATA ARRAY] 
if "`mata'"=="" & "`array'"==""{
    return local ldot=","
}

end


capture program drop parsespmat1
program define parsespmat1, rclass
syntax namelist(name=wnames),aname(name) [MATA ARRAY] 
local nw: word count `wnames'

local i=1
if "`mata'"=="" & "`array'"==""{
	mata: `aname' = asarray_create("real")
	foreach w in `wnames'{
		tempname w`i' w`i'_id
		spmatrix matafromsp `w`i'' `w`i'_id' = `w'
		local matanames `matanames' `w`i''
		mata: asarray(`aname',`i',`w`i'')
		local i=`i'+1
	}
	cap mata mata drop `matanames'
	
}
else if "`mata'"!=""{
	mata: `aname' = asarray_create("real")
	local matanames `wnames'
	local i=1
	foreach w in `matanames'{
		mata: asarray(`aname',`i',`w')
		local i=`i'+1
	}

}
else{
	mata: _temparray = asarray_create("real")
	mata: keys = asarray_keys(`wnames')
    mata: keys = sort(keys,1) // sort w in time order
	mata: st_local("keytypes",eltype(keys))
	if ("`keytypes'"!="real"){
		di as error "keys in array `wnames' is not real"
		exit 198
	}
    mata: st_numscalar("r(nw)",length(keys))
	local nw = r(nw)
    forv j=1/`nw'{
		mata: asarray(_temparray,`j',asarray(`wnames',keys[`j']))
	}
	mata: `aname' = asarray_create("real")
	forv j=1/`nw'{
		mata: asarray(`aname',`j',asarray(_temparray,`j'))
	}
	cap mata mata drop _temparray

 }

return scalar nw=`nw'

end 


//////////////////////////////////////
* This command is browed from Morad Zekhnini's nwxtregress package
capture program drop issorted
program define issorted
	syntax	varlist 
	
	local sorted : dis "`:sortedby'"
	if "`sorted'" != "`varlist'" {
	    noi disp `"sort data by `varlist'"'
		noi disp "make sure that each spmatrix is the same order" _n
	    sort `varlist'
	}

end

// -------------------------------------------------------------------------------------------------

///////
capture program drop data2mata
program define data2mata,rclass
version 16

gettoken mname 0: 0, p(=)
gettoken eq 0:0 

if usubinstr("`0'"," ","",.)==""{
	mata: `mname' = 1
	return scalar ncol = 1
	exit
}

syntax varlist(min=1) [if] [in], [NOConstant] 
marksample touse

mata: `mname' = st_data(., `"`varlist'"', "`touse'")
if "`noconstant'"==""{
	mata: `mname' = `mname', J(rows(`mname'),1,1)
}
mata: st_numscalar("r(nc)",cols(`mname'))
return scalar ncol = r(nc)


end

// -------------------------------------------------------------------------------------------------

cap mata mata drop spte()
cap mata mata drop sphut()
mata:
/*
void function spte(transmorphic matrix w,
                    real scalar rho,
				    real scalar T,
					real colvector te_order,
					real colvector htildeut)
{
	external real colvector _pan_tvar
	external real colvector _pan_tvar
	info = panelsetup(_pan_tvar,1)
	if (rho==.){
		htildeut = ln(te_order)
	}
	else{
		keys = asarray_keys(w)
		htildeut = J(0,1,.)
		if (length(keys)==1){
			wt = asarray(w,1)
			irhow = matinv(I(rows(wt))-rho*wt)
		}

		for(t=1;t<=T;t++){

			if(length(keys)>1){
				wt = asarray(w,t)
				irhow = matinv(I(rows(wt))-rho*wt)
			}
			htildeut = htildeut \ irhow*ln(panelsubmatrix(te_order,t,info))
		}
		te_order = exp(htildeut)
	}

}
*/

real colvector spte(real matrix xx, 
					real colvector bx,
					real matrix vv, 
					real colvector bv,
					real matrix uu, 
					real colvector bu,
					real scalar rho, 
					real scalar gamma, 
					real scalar tau,
					real colvector y ,
					transmorphic matrix wy,
					transmorphic matrix wv,
				    transmorphic matrix wu)
			                              
{
	
external real scalar _cost
external real colvector _pan_tvar
//external real matrix info
info = panelsetup(_pan_tvar,1)
//external transmorphic matrix wy_ina
//external transmorphic matrix wv_ina
//external transmorphic matrix wu_ina
nt = panelstats(info)[1]
mu = 0
es = _cost*(y - xx*bx )
if (vv==1){
	sigv2 = exp(bv)
}
else{
	sigv2 = exp(vv*bv)
}
hb = exp(0.5*uu*bu)
//s2_u	= exp(b[nx+nz+2])	// scalar
s2_u=1

if(wy!=.){
	wykeys = asarray_keys(wy)
	lenwykeys = length(wykeys)
}
else{
	lenwykeys = 0
}
if(wv!=.){
	wvkeys = asarray_keys(wv)
	lenwvkeys = length(wvkeys)
}
else{
	lenwvkeys = 0
}
if(wu!=.){
	wukeys = asarray_keys(wu)
	lenwukeys = length(wukeys)
}
else{
	lenwukeys = 0
}
lndetPi =.
invPi =.
if (lenwvkeys==1){
	spwi = asarray(wv,1)
	Mr = I(rows(spwi))-rho*spwi
	iMr = matinv(Mr)
	pifun(sigv2,Mr,iMr,lndetPi,invPi)
}
else{
	Mr = 1
	iMr = 1
	pifun(sigv2,Mr,iMr,lndetPi,invPi)
}
if (lenwukeys==1){
	spwi = asarray(wu,1)
	Mtau = matinv(I(rows(spwi))-tau*spwi)
}
else{
	Mtau = 1
}

if(lenwykeys==1){
	wyi = asarray(wy,1)
	Mrho = matinv(I(rows(wyi))-rho*wyi)
}
else{
	Mrho = 1
}

Eie = J(0,1,.)
for(i=1;i<=nt;i++){
	N = info[i,2]-info[i,1]+1
	if(lenwukeys>1){
		spwi = extrpoi(asarray(wu,i))
		Mtau = matinv(I(N)-tau*spwi)          
	}
	if(lenwvkeys>1){
		spwi = extrpoi(asarray(wv,i))
		Mr = matinv(I(N)-gamma*spwi)
		iMr = matinv(Mr)
		pifun(sigv2,Mr,iMr,lndetPi,invPi)
	}

	tvpifun(sigv2,Mr,iMr,lndetPi,invPi,i,info)

	hbi	= Mtau*panelsubmatrix(hb, i, info)
	esi = panelsubmatrix(es, i, info)
	hinvPi = quadcross(hbi,invPi)
	hinvPihi = quadcross(hinvPi',hbi)
	esinvPi = quadcross(esi,invPi)
	s2_s = 1/(hinvPihi+1/s2_u)
	us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
	htildeuts = hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))	
   if(lenwykeys>1){
		wyi = asarray(wy,i)
		Mrho = matinv(I(N)-rho*wyi)
	}
	htildeuts = Mrho*htildeuts
	Eie = Eie \ htildeuts
}

return(exp(-Eie))
	
}


// compute the estimate h(z*dela)*E(u_t^*|e_t*) for h(z*dela)*u_t^*
real colvector sphut(real matrix xx, 
					real colvector bx,
					real matrix vv, 
					real colvector bv,
					real matrix uu, 
					real colvector bu,
					real scalar rho, 
					real scalar gamma, 
					real scalar tau,
					real colvector y ,
					transmorphic matrix wv,
					transmorphic matrix wu)
						  
{

external real scalar _cost
external real colvector _pan_tvar
//external real matrix info
info = panelsetup(_pan_tvar,1)
//external transmorphic matrix wy_ina
//external transmorphic matrix wv_ina
//external transmorphic matrix wu_ina
nt = panelstats(info)[1]
mu = 0
es = _cost*(y - xx*bx )
if (vv==1){
	sigv2 = exp(bv)
}
else{
	sigv2 = exp(vv*bv)
}
hb = exp(0.5*uu*bu)
//s2_u	= exp(b[nx+nz+2])	// scalar
s2_u=1

if(wv!=.){
	wvkeys = asarray_keys(wv)
	lenwvkeys = length(wvkeys)
}
else{
	lenwvkeys = 0
}
if(wu!=.){
	wukeys = asarray_keys(wu)
	lenwukeys = length(wukeys)
}
else{
	lenwukeys = 0
}
lndetPi =.
invPi =.
if (lenwvkeys==1){
	spwi = asarray(wv,1)
	Mr = I(rows(spwi))-rho*spwi
	iMr = matinv(Mr)
	pifun(sigv2,Mr,iMr,lndetPi,invPi)
}
else{
	Mr = 1
	iMr = 1
	pifun(sigv2,Mr,iMr,lndetPi,invPi)
}
if (lenwukeys==1){
	spwi = asarray(wu,1)
	Mtau = matinv(I(rows(spwi))-tau*spwi)
}
else{
	Mtau = 1
}
Eie = J(0,1,.)
for(i=1;i<=nt;i++){
	N = info[i,2]-info[i,1]+1
	if(lenwukeys>1){
		spwi = extrpoi(asarray(wu,i))
		Mtau = matinv(I(N)-tau*spwi)          
	}
	if(lenwvkeys>1){
		spwi = extrpoi(asarray(wv,i))
		Mr = matinv(I(N)-gamma*spwi)
		iMr = matinv(Mr)
		pifun(sigv2,Mr,iMr,lndetPi,invPi)
	}
	esi = panelsubmatrix(es, i, info)
	tvpifun(sigv2,Mr,iMr,lndetPi,invPi,i,info)
	hbi0 = panelsubmatrix(hb, i, info)
	hbi	= Mtau*hbi0
	hinvPi = quadcross(hbi,invPi)
	hinvPihi = quadcross(hinvPi',hbi)
	esinvPi = quadcross(esi,invPi)
	s2_s = 1/(hinvPihi+1/s2_u)
	us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
	htildeuts = hbi0*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))	
	Eie = Eie \ htildeuts
}

return(Eie)

}

end
