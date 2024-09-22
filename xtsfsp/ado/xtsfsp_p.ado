*! 2024-07-04

cap program drop xtsfsp_p
cap program define xtsfsp_p, rclass sortpreserve
version 16
syntax newvarname, [xb Residuals uts u su te ste sigv hs]
local nopts: word count `xb' `residuals' `uts' `u' `te'  `su'  `ste' `sigmav' `hs' 
    if `nopts' >1 {
        display "{err}only one statistic may be specified"
        exit 498
    }
	if `nopts'==0 local xb xb

// extract ereturn information from last xtsfsp estimation
local cmd `e(cmd)'
if ustrpos("`cmd'", "xtsfsp") == 0 {
    di as err "last estimate of xtsfsp not found"
    exit
}
local depvar `e(depvar)'
local rho = e(rho)
local gamma = e(gamma)
local tau = e(tau)
local wymat  `e(wy)'
local wvmat  `e(wv)'
local wumat  `e(wu)'
local hasgenwvars `e(hasgenwvars)'
local ivar `e(ivar)'
local tvar `e(tvar)'
local veq  `e(veq)'

scalar rymin = e(rymin)
scalar rymax = e(rymax)
scalar rumin = e(rumin)
scalar rumax = e(rumax)
scalar rvmin = e(rvmin)
scalar rvmax = e(rvmax)

if(usubinstr("`hasgenwvars'"," ","",.)==""){
  di as error "genwvars must be specified in xtsfsp estimate for using the postestimation command"
  exit
}

tempvar yhat vhet uhet ehat sigmav2 ht

qui _predict double `yhat', xb eq(#1)
if (`"`rho'"'!="" & `"`rho'"'!="." & `"`rho'"'!="0"){
	qui replace `yhat' = `yhat' + `rho'*W_`depvar'
}  
if "`xb'"!=""{
    qui gen double `varlist' = `yhat'
    label var `varlist' "prediction of `depvar'"
    exit
}  
qui gen double `ehat' = `depvar' - `yhat'
if "`residuals'"!="" {
    qui gen double `varlist' = `ehat'
    label var `varlist' "residuals"
    exit
}
qui _predict double `vhet', xb eq(#2)
qui gen double `sigmav2' = exp(`vhet')

if "`sigmav'"!=""{
    qui gen double `varlist' = sqrt(exp(`sigmav2'))
    label var `varlist' "sigma_v"
    exit
}
qui _predict double `uhet', xb eq(#3)
qui gen double `ht' = exp(0.5*`uhet')
if "`hs'"!=""{
    qui gen double `varlist' = `ht'
    label var `varlist' "estimate of scaling function"
    exit
}
if ("`u'"!="" | "`su'"!=""){
   confirm new var `varlist'_dir
   confirm new var `varlist'_indir
}

tempvar time2 order0
qui egen `time2' = group(`tvar')
sort `time2' `ivar'

preserve
qui gen `order0' = _n
cap keep if __e_sample__
mata: _order_0 = st_data(., "`order0'")
mata:`ehat'= st_data(.,"`ehat'")
//mata: `sigmav2' = st_data(.,"`sigmav2'")
mata: `ht' = st_data(.,"`ht'")
if(usubinstr("`veq'"," ","",.)==""){ // only a constant 
    local sigv2 = `sigmav2'[1]
    mata: `sigmav2' = `sigv2'
}
else{
    mata: `sigmav2' = st_data(., "`sigmav2'") 
}

if "`uts'"!=""{
    mata:_u_order_= sputs(`ehat', `sigmav2', `ht', `rho', `gamma', `tau', _cost, parray(`wvmat'), parray(`wumat'))
}

if "`u'"!=""| "`te'"!=""{
    mata:_u_order_= sputl(`ehat', `sigmav2', `ht', `rho', `gamma', `tau', _cost, parray(`wvmat'), parray(`wumat'), _dirsu=., _indirsu=.)
}
if "`su'"!=""| "`ste'"!=""{
    mata:_u_order_= spsu(`ehat', `sigmav2', `ht', `rho', `gamma', `tau', _cost, parray(`wymat'), parray(`wvmat'), parray(`wumat'), _dirsu=., _indirsu=.)
}

restore
/* put back estimates into origin data*/
if "`uts'"!=""{
    qui gen double `varlist' = .
    label var `varlist' "u*_t estimate"
    mata: getdatafmata(_u_order_,_order_0,"`varlist'")
    cap mata mata drop  _u_order_
    exit
}
if "`u'"!=""{
    qui gen double `varlist' = .
    label var `varlist' "inefficiency estimate E(u|e)"
    mata: getdatafmata(_u_order_,_order_0,"`varlist'")
    cap mata mata drop  _u_order_
    qui gen double `varlist'_dir = .
    label var `varlist'_dir "inefficiency estimate E(u|e): direct"
    mata: getdatafmata(_dirsu, _order_0,"`varlist'_dir")
    cap mata mata drop  _dirsu
    qui gen double `varlist'_indir = .
    label var `varlist'_indir "inefficiency estimate E(u|e): indirect"
    mata: getdatafmata(_indirsu, _order_0,"`varlist'_indir")
    cap mata mata drop  _indirsu
    exit
}
if "`te'"!=""{
    qui gen double `varlist' = .
    label var `varlist' "efficiency estimate exp(-E(u|e))"
    mata: _u_order_= exp(-_u_order_)
    mata: getdatafmata(_u_order_,_order_0,"`varlist'")
    cap mata mata drop  _u_order_
    exit
}

if "`su'"!=""{
    qui gen double `varlist' = .
    label var `varlist' "spatial inefficiency estimate [I-rh0*W]^{-1}E(u|e)"
    mata: getdatafmata(_u_order_,_order_0,"`varlist'")
    cap mata mata drop  _u_order_
    qui gen double `varlist'_dir = .
    label var `varlist'_dir "spatial inefficiency estimate [I-rh0*W]^{-1}E(u|e): direct"
    mata: getdatafmata(_dirsu, _order_0,"`varlist'_dir")
    cap mata mata drop  _dirsu
    qui gen double `varlist'_indir = .
    label var `varlist'_indir "spatial inefficiency estimate [I-rh0*W]^{-1}E(u|e): indirect"
    mata: getdatafmata(_indirsu, _order_0,"`varlist'_indir")
    cap mata mata drop  _indirsu
    exit
}
if "`ste'"!=""{
    qui gen double `varlist' = .
    label var `varlist' "spatial efficiency estimate exp(-[I-rh0*W]^{-1}E(u|e))"
    mata: _u_order_= exp(-_u_order_)
    mata: getdatafmata(_u_order_,_order_0,"`varlist'")
    cap mata mata drop  _u_order_
    exit
}

end

cap mata mata drop parray()
cap mata mata drop spsu()
cap mata mata drop sputl()
cap mata mata drop sputs()
mata:
real colvector spsu(real colvector ehat, 
					real colvector sigv2,
					real colvector hb,
					real scalar rho, 
					real scalar gamma, 
					real scalar tau,
                    real scalar _cost,
					transmorphic matrix wy,
					transmorphic matrix wv,
				    transmorphic matrix wu,
                    real colvector diru,
                    real colvector idiru)
			                              
{
external real colvector _pan_tvar
//external real matrix info
info = panelsetup(_pan_tvar,1)
//external transmorphic matrix wy_ina
//external transmorphic matrix wv_ina
//external transmorphic matrix wu_ina
nt = panelstats(info)[1]
mu = 0
es = _cost*ehat
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
diru =J(0,1,.)
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
	else{
		invPi = I(N)*invPi
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
    hbidu = diag(diagonal(Mrho*Mtau))*panelsubmatrix(hb, i, info)
	Eie = Eie \ htildeuts
    diru = diru \ hbidu*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))
}

idiru = Eie - diru

return(Eie)
	
}


real colvector sputl(real colvector ehat, 
					real colvector sigv2,
					real colvector hb,
					real scalar rho, 
					real scalar gamma, 
					real scalar tau,
                    real scalar _cost,
					transmorphic matrix wv,
				    transmorphic matrix wu,
                    real colvector diru,
                    real colvector idiru)
			                              
{
external real colvector _pan_tvar
//external real matrix info
info = panelsetup(_pan_tvar,1)
//external transmorphic matrix wy_ina
//external transmorphic matrix wv_ina
//external transmorphic matrix wu_ina
nt = panelstats(info)[1]
mu = 0
es = _cost*ehat
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
	//invPi
}
if (lenwukeys==1){
	spwi = asarray(wu,1)
	Mtau = matinv(I(rows(spwi))-tau*spwi)
}
else{
	Mtau = 1
}


Eie = J(0,1,.)
diru =J(0,1,.)
//indiru = J(0,1,.)
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
	else{
		invPi = I(N)*invPi
	}

	tvpifun(sigv2,Mr,iMr,lndetPi,invPi,i,info)

	hbi	= Mtau*panelsubmatrix(hb, i, info)

    hbidu = diag(diagonal(Mtau))*panelsubmatrix(hb, i, info)
	esi = panelsubmatrix(es, i, info)
	//rows(hbi),rows(invPi),rows(esi)
	hinvPi = quadcross(hbi,invPi)
	hinvPihi = quadcross(hinvPi',hbi)
	esinvPi = quadcross(esi,invPi)
	
	s2_s = 1/(hinvPihi+1/s2_u)
	us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
	htildeuts = hbi*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))	
	Eie = Eie \ htildeuts
    diru = diru \ hbidu*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))	
}
    idiru = Eie - diru

return(Eie)
	
}

real colvector sputs(real colvector ehat, 
					real colvector sigv2,
					real colvector hb,
					real scalar rho, 
					real scalar gamma, 
					real scalar tau,
                    real scalar _cost,
					transmorphic matrix wv,
				    transmorphic matrix wu)
			                              
{
external real colvector _pan_tvar
//external real matrix info
info = panelsetup(_pan_tvar,1)
//external transmorphic matrix wy_ina
//external transmorphic matrix wv_ina
//external transmorphic matrix wu_ina
nt = panelstats(info)[1]
mu = 0
es = _cost*ehat
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
	else{
		invPi = I(N)*invPi
	}

	tvpifun(sigv2,Mr,iMr,lndetPi,invPi,i,info)

	hbi	=Mtau* panelsubmatrix(hb, i, info)
	esi = panelsubmatrix(es, i, info)
	hinvPi = quadcross(hbi,invPi)
	hinvPihi = quadcross(hinvPi',hbi)
	esinvPi = quadcross(esi,invPi)
	s2_s = 1/(hinvPihi+1/s2_u)
	us = (mu/s2_u - quadcross(esinvPi',hbi))*s2_s
	htildeuts = J(N,1,1)*(us+sqrt(s2_s)*normalden(us/sqrt(s2_s))/normal(us/sqrt(s2_s)))	
	Eie = Eie \ htildeuts
}

return(Eie)
	
}


function parray(transmorphic matrix w)
{
		if(eltype(w)=="pointer"){
			return(*w)
		}
		else{
			return(w)
		}
}

end
