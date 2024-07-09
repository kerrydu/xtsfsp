*! 2024-07-03
*! Monday, July 1, 2024 at 13:22:19, by Kerry Du
*! total, direct and indirect effects for xtsfsp
capture program drop xtsfsp_margins
program define xtsfsp_margins, rclass

version 16
syntax varlist, [reps(real 200) seed(real 123) NODOTS NORMalize FIXUts]

    if "`nodots'"!="" local qui qui
    // default: compute absolute ME
    if "`normalize'"=="" local absolute absolute
    local marginvars `varlist'
    ////////////////////////////////////////
    * extract informtaion from estimates
    local cmd  `e(cmd)'
    local cmdline `e(cmdline)'
    local yvar `e(depvar)'
    local wymat `e(wy)'
    local wxmat `e(wx)'
    local wumat `e(wu)'
    local wvmat `e(wv)'
    local wxvars `e(wxvars)'
    local T = e(T)
    local rho = e(rho)
    local gamma =e(gamma)
    local tau = e(tau)
    local xeq  `e(xeq)'
    local veq `e(veq)'
    local ueq `e(ueq)'
    local hasgenwvars `e(hasgenwvars)'

    // count parameters in the first 3 eqs
    if `"`xeq'"'=="" | `"`xeq'"'=="."{ // only has a constant
        local nx = 1
    }
    else{
        countnumv `xeq'
        local nx = r(nv)
    }

    if `"`veq'"'=="" | `"`veq'"'=="."{ // only has a constant
        local nv=1
    }
    else{
        countnumv `veq'
        local nv = r(nv)
    }
    if `"`ueq'"'=="" | `"`ueq'"'=="."{ // only has a constant
        local nz=1
    }
    else{
        countnumv `ueq'
        local nz = r(nv)
    }

    if "`hasgenwvars'"==""{
        di as err "genwvars MUST be specified in the last estimates"
        exit
    }

    scalar rymin = e(rymin)
    scalar rymax = e(rymax)
    scalar rumin = e(rumin)
    scalar rumax = e(rumax)
    scalar rvmin = e(rvmin)
    scalar rvmax = e(rvmax)
    tempname bml covml
    mat `bml'  = e(b)
    mat `covml' = e(V)
    mata: VCV = st_matrix("`covml'")
    mata: VCV = cholesky(VCV)
    mata: bml = st_matrix("`bml'")
	local varnames: colnames `bml'
    local nparas : word count `varnames'
	mata: varnames=tokens(st_local("varnames"))
    gettoken cmd0 cmdline : cmdline
    cmdlineparse `cmdline' // parse the cmdline of last estimate
    local wy0 `r(wy0)' // original wy information
    local wu0 `r(wu0)' // original wu information
    local wv0 `r(wv0)' // original wv information
    local taucon = .
    local rhocon = .
    local rhotauflag = 0
    local tau_index = 0
    local rho_index = 0
    local gamma_index =0
    local xvars = r(xvars)
    local zvars = r(zvars)
    if `"`wu0'"'!= "" {
        local tau_index = `nparas' // para = bx, b_vhet, b_uhet, rho, theta,
    }
    else{
        local tau = 0
        local wumat wuina 
        mata: wuina  = asarray_create("real")
        mata: asarray(wuina,1,1) 
    }
    if `"`wy0'"'!= "" {
        local rho_index = `nparas' // para = bx, b_vhet, b_uhet, rho, theta, tau
        if `"`wu0'"'!= "" {
            local rho_index = `rho_index' -1
        }
        if `"`wv0'"'!= "" {
            local rho_index = `rho_index' -1
        }
    }
    else{
        local rho = 0
        local wymat wyina 
        mata: wyina  = asarray_create("real")
        mata: asarray(wyina,1,1) 
    }
    if `"`wv0'"'!=""{
        local gamma_index = `nparas' // para = bx, b_vhet, b_uhet, rho, theta, tau
        if `"`wu0'"'!= "" {
            local gamma_index = `gamma_index' -1
        }
    }
    else{
        local gamma = 0
        local wvmat wvina 
        mata: wvina  = asarray_create("real")
        mata: asarray(wvina,1,1) 
    }

    **********************
    if ustrpos("`cmd'", "xtsfsp") == 0 {
        di as err "xtsfsp_margins: last estimate of xtsfsp not found"
        exit
    }

	local allvars `xvars' `zvars' 
	local ccom: list marginvars - allvars 
	if `"`ccom'"'!=""{
		di as error `"variables {`ccom'} not specified in the last estimates"'
        di as red "Check variables specified in frontier and uhet() of the last estimates"
		exit 
	}
	local ccomx: list marginvars & xvars
	local ccomu: list marginvars & zvars     

    confirm  integer number `reps'
    confirm integer number `seed'
    if `reps' < 0{
        di as err "reps() should be a positive integer"
        exit
    }

   preserve
   //!!!sort time id!!!
   cap confirm var  __e_sample__ 
   if _rc ==0{
     keep if __e_sample__==1
   }
   _xt, trequired 
   local id=r(ivar)
   local time=r(tvar)
   tempvar time2
   qui egen `time2' = group(`time')
   sort `time2' `id'

   mata: _pan_tvar = st_data(., "`time2'")
//    mata: hzt = st_data(., "`hzt'") 
   mata: transinvw(`wymat',`T',`rho',irhow=.)

if `"`ccomx'"'!=""{
    //di "Marginal effects for variables in the frontier: "
    local nmx: word count `ccomx'
    mata: sdxt = J(`nmx',0,.)
    mata: sdxd = J(`nmx',0,.)
    mata: sdxi = J(`nmx',0,.)

    if ("`wxmat'"==""| "`wxmat'"=="."){
        mata: wxina  = asarray_create("real")
        mata: wxina = asarray(wxina,1,0)
        local wxmat wxina
    }

    mata:totexmat=x_mc(`"`ccomx'"',varnames,bml',`T',irhow,`wxmat',direxmat=.,indirexmat=.)
}


if `"`ccomu'"'!=""{
   // di "Marginal effects for variables in the inefficiency function: "
   local nmx: word count `ccomu'
    mata: sdzt = J(`nmx',0,.)
    mata: sdzd = J(`nmx',0,.)
    mata: sdzi = J(`nmx',0,.)
    mata: totezmat = J(0,1,.)
    mata: direzmat = J(0,1,.)
    mata: indirezmat = J(0,1,.)
    mata: transinvw(`wumat',`T',`tau',itauw=.)
    tempvar hzt uts zzz xxx yyy vvv 
    data2mata `zzz' = `ueq'
    xtsfsp_p `hzt', hs 
    if "`absolute'"!=""{
        xtsfsp_p `uts', uts 
        qui replace `hzt' = `hzt'*`uts' 
        mata: uts = st_data(., "`uts'")
        mata: `yyy' = st_data(.,"`yvar'")
        data2mata `xxx' = `xeq'
        if(`"`veq'"'=="") local veq = 1
        else data2mata `vvv' = `veq'
        if(`rho'==0 | `rho'==.){
            local wyvar = 1
        }
        else{
            tempname wyvar 
            mata: `wyvar' = st_data(.,"W_`yvar'")
        }
    }
    mata: hzust = st_data(., "`hzt'")
    
    local nxv =`nx' + `nv'
    mata: bz=bml[(`nxv'+1)..(`nxv'+`nz')]
    mata:totezmat= z_mc(`"`ccomu'"',varnames[`=`nxv'+1'..length(varnames)],bz',`T',hzust,irhow,itauw,direzmat,indirezmat)


}

** mc to compute standard errors
set seed `seed'

if `reps'>0 di _n "Monte Carlo simulation: " 
mata: bml0 = bml'
forv b=1/`reps'{
    `qui' displaydot , dot(`b')
    mata: bml =bml0+ bmc(VCV)
    //mata: bml
    mata: rho_gamma_tau(bml,`rho_index',`gamma_index',`tau_index',rho=.,gamma=.,tau=.,1)
    mata: transinvw(`wymat',`T',rho,irhow=.)
    if `"`ccomx'"'!=""{
        mata: totexb = x_mc(`"`ccomx'"',varnames,bml,`T',irhow,`wxmat',direxb=.,indirexb=.)
        mata: sdxt = sdxt, totexb
        mata: sdxd = sdxd, direxb
        mata: sdxi = sdxi, indirexb
    }
    
    if `"`ccomu'"'!=""{
        if `tau'==0 mata: itau = .
        else mata: transinvw(`wumat',`T',tau,itauw=.)
        mata: bz=bml[(`nxv'+1)..(`nxv'+`nz')]
        if "`absolute'"==""{
            mata: hzust = exp(0.5*`zzz'*bz)
        }
        else if("`fixuts'"==""){
            if `gamma'==0 mata: igammaw = .
            else mata: transinvw(`wvmat',`T',gamma,igammaw=.)
            mata: hzust = rst_hzust(`yyy',`wyvar',`xxx',`vvv',`zzz',bml,rho,igammaw,itauw)
        }
        else{
            mata: hzust = exp(0.5*`zzz'*bz):*uts
        }

        mata: bz=bml[(`nxv'+1)..(`nxv'+`nz')]
        mata: totezb= z_mc(`"`ccomu'"',varnames[`=`nxv'+1'..length(varnames)],bz,`T',hzust,irhow,itauw,direzb=.,indirezb=.)
        mata: sdzt = sdzt, totezb
        mata: sdzd = sdzd, direzb
        mata: sdzi = sdzi, indirezb
   }
    
}

 if `"`ccomx'"'!=""{
    //mata: sdxt',sdxd',sdxi'
    mata: sdxt = rowsd(sdxt)
    mata: sdxd = rowsd(sdxd)
    mata: sdxi = rowsd(sdxi)
    mata: totexmat =totexmat,sdxt, (totexmat:/sdxt), 2*(1:-normal(abs(totexmat:/sdxt)))
    mata: direxmat =direxmat, sdxd,(direxmat:/sdxd), 2*(1:-normal(abs(direxmat:/sdxd)))
    mata: indirexmat = indirexmat,sdxi, (indirexmat:/sdxi), 2*(1:-normal(abs(indirexmat:/sdxi)))   
    di _n "Marginal effects for variables in the frontier: "
    mata: st_matrix("totxeff",totexmat)
    local rnames
    foreach v in `ccomx'{
        local rnames `rnames' `"`v'"'
    }
    mat rownames totxeff = `rnames'
    mat colnames totxeff = "Coeff" "se" "z" "P"			//mata: _b_ml = st_matrix("`bml'")

    local r = rowsof(totxeff)-1
    local rf "--"
    forvalues i=1/`r' {
        local rf "`rf'&"
    }
    local rf "`rf'-"
    local cf "&  %10s | %12.4f & %12.4f &  %12.4f  & %12.4f &"
    dis _n in gr "Total marginal effects:"
    matlist totxeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")			

    mata: st_matrix("direxeff",direxmat)
    mat rownames direxeff = `rnames'
    mat colnames direxeff = "Coeff" "se" "z" "P"			
    dis _n in gr "Direct marginal effect:"
    matlist direxeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
    mata: st_matrix("indirexeff",indirexmat)
    mat rownames indirexeff = `rnames'
    mat colnames indirexeff = "Coeff" "se" "z" "P"		
    dis _n in gr "Indirect marginal effect:"
    matlist indirexeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")

    return matrix totalxmargins = totxeff
    return matrix directxmargins = direxeff
    return matrix indirectxmargins = indirexeff

 }

 if `"`ccomu'"'!=""{
   // mata: sdzt', sdzd', sdzi'
    mata: sdzt = rowsd(sdzt)
    mata: sdzd = rowsd(sdzd)
    mata: sdzi = rowsd(sdzi)
    mata: totezmat =totezmat,sdzt, (totezmat:/sdzt), 2*(1:-normal(abs(totezmat:/sdzt)))
    mata: direzmat =direzmat,sdzd, (direzmat:/sdzd), 2*(1:-normal(abs(direzmat:/sdzd)))
    mata: indirezmat = indirezmat,sdzi, (indirezmat:/sdzi), 2*(1:-normal(abs(indirezmat:/sdzi)))
    di _n "Marginal effects for variables in the inefficiency function: "
    mata: st_matrix("totzeff",totezmat)
    local rnames
    foreach v in `ccomu'{
        local rnames `rnames' `"`v'"'
    }
    mat rownames totzeff = `rnames'
    mat colnames totzeff = "Coeff" "se" "z" "P"			//mata: _b_ml = st_matrix("`bml'")

    local r = rowsof(totzeff)-1
    local rf "--"
    forvalues i=1/`r' {
        local rf "`rf'&"
    }
    local rf "`rf'-"
    local cf "&  %10s | %12.4f & %12.4f &  %12.4f  & %12.4f &"
    dis _n in gr "Total marginal effects:"
    matlist totzeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
                
    mata: st_matrix("direzeff",direzmat)
    mat rownames direzeff = `rnames'
    mat colnames direzeff = "Coeff" "se" "z" "P"			
    dis _n in gr "Direct marginal effect:"
    matlist direzeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")

    mata: st_matrix("indirezeff",indirezmat)
    mat rownames indirezeff = `rnames'
    mat colnames indirezeff = "Coeff" "se" "z" "P"			
    dis _n in gr "Indirect marginal effect:"
    matlist indirezeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
    return matrix totalzmargins = totzeff
    return matrix directzmargins = direzeff
    return matrix indirectzmargins = indirezeff


 }


di  _n "Note: Standard errors are computed using the Monte Carlo method."
di     "      P value are computed using normal distribution."
restore

end



*************************************************************************************

cap program drop displaydot
 program define displaydot
 version 14
 syntax, dot(integer)
    if mod(`dot',100) == 0 {
        di _c `dot' 
        di
    }
    else if mod(`dot',5) == 0{
        disp _c "."
    }
  end

****************************************************************************************
cap program drop cmdlineparse
program define cmdlineparse,rclass
version 16

syntax varlist,Uhet(string) [wy(string) wx(string) wu(string) wv(string) *] 

gettoken yvar xvar: varlist
if (`"`xvar'"'!="" | "`nocostant'"!="")  {
    countnumv `xvar',`noconstant'
    local nx = r(nv)
}
else{
    local nx = 1
}
if (`"`uhet'"'!=""){
    countnumv `uhet'
    local nu = r(nv)
}
else{
    local nu =1
}
if (`"`vhet'"'!=""){
    countnumv `vhet'
    local nv = r(nv)
}
else{
    local nv = 1
}

gettoken zvars : uhet, p(,)
gettoken vvars :vhet, p(,) //v

return scalar nx = `nx'
return scalar nu = `nu'
return scalar nv = `nv'
return local zvars  `zvars'
return local vvars  `vvars'
return local xvars  `xvar'
return local yvar  `yvar'
return local wy0 `wy'
return local wu0 `wu'
return local wv0 `wv'

end

cap program drop countnumv
program define countnumv, rclass
	version 16
	syntax [varlist],[NOCONstant]
    if(`"`varlist'"'=="" & "`noconstant'"!=""){
        di as err "No variables are specified"
        exit
    }
    if(`"`varlist'"'=="" & "`noconstant'"==""){
        return scalar nv = 1
        exit
    }    
	local nv: word count `varlist'
	if "`noconstant'"==""{
		local nv = `nv'+1
	}
	return scalar nv = `nv'
end

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



*****************************************************************************

/*!
 * @author Kerry Du
 * @date 2024-07-01
 * @version 1.0
 * @brief This file contains the mata functions for computing the total, direct and indirect effects for xtsfsp
 * @details The total, direct and indirect effects are computed for the frontier and inefficiency function
 * @see xtsfsp_margins_delta, xtsfsp_margins_mc
 */
 
 cap mata mata drop rst_hzust()
 cap mata mata drop x_mc()
 cap mata mata drop z_mc()
 cap mata mata drop marginx_mc()
 cap mata mata drop marginz_mc() 
 cap mata mata drop transinvw()
 cap mata mata drop extractposx()
 cap mata mata drop rowsd()
 cap mata mata drop bmc()
 cap mata mata drop rho_tau()
 cap mata mata drop rho_gamma_tau()

 
mata:
// estimate h(z)*E(u_t^*|e_t)
real colvector rst_hzust(real colvector yy,
                         real colvector wyy,
                         real matrix xx,
                         real matrix vv,
                         real matrix zz,
                         real colvector bb,
                         real scalar rho,
                         transmorphic matrix igammaw,
                         transmorphic matrix itauw)
						  
{
external real scalar _cost
external real colvector _pan_tvar
bx = bb[1::cols(xx)]
bv = bb[(cols(xx)+1)::(cols(xx)+cols(vv))] 
bz = bb[(cols(xx)+cols(vv)+1)::(cols(xx)+cols(vv)+cols(zz))]   

//external real matrix info
info = panelsetup(_pan_tvar,1)
//external transmorphic matrix wy_ina
//external transmorphic matrix wv_ina
//external transmorphic matrix wu_ina
nt = panelstats(info)[1]
mu = 0
es = _cost*((yy - xx*bx):-rho*wyy)
if (vv==1){
	sigv2 = exp(bv)
}
else{
	sigv2 = exp(vv*bv)
}
hb = exp(0.5*zz*bz)
//s2_u	= exp(b[nx+nz+2])	// scalar
s2_u=1

lndetPi =.
invPi =.


if(igammaw!=.){
    wvkeys = asarray_keys(igammaw)
    lenwvkeys = length(wvkeys)
    if(length(wvkeys)==1){
        Mr = asarray(igammaw,wvkeys[1])
        iMr = matinv(Mr)
        pifun(sigv2,Mr,iMr,lndetPi,invPi)
    }
}
else{
    lenwvkeys = 0
    Mr = 1
    iMr = 1
    pifun(sigv2,Mr,iMr,lndetPi,invPi)
}

if(itauw!=.){
    wukeys = asarray_keys(itauw)
    lenwukeys = length(wukeys)
    if(length(wukeys)==1){
        Mtau = asarray(itauw,wukeys[1])
    }
}
else{
    lenwukeys = 0
    Mtau = 1
}


Eie = J(0,1,.)

for(i=1;i<=nt;i++){
	//N = info[i,2]-info[i,1]+1
	if(lenwukeys>1){
		Mtau = extrpoi(asarray(itauw,i))         
	}
	if(lenwvkeys>1){
		Mr = extrpoi(asarray(igammaw,i))
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

// find postion of of coefficienct of x and w_x in the parameter vector
real vector extractposx(string rowvector varnames, string scalar var)
{ 
	 k = select(1..length(varnames),varnames:==var)
	 k1 = k[1]
	 k = select(1..length(varnames),varnames:==("W_"+var))
	if(length(k) >0){
		k2 = k[1]	
		k1 = k1, k2
	}
	return(k1)
}

// tote=marginx(`rho',bi,VCV[ii,ii],`T',w_ina,`wxmat',dire=.,indire=.,)
real scalar marginx_mc( real rowvector b, 
                        real scalar T, 
                        transmorphic matrix iwy, 
                        transmorphic matrix wx, 
                        real scalar dire, 
                        real scalar idire)
{
    keys = asarray_keys(iwy)
    if(length(keys)==1){
        irhow = asarray(iwy,keys[1])
        if (length(b)==1){
            totale =sum(irhow*b)/rows(irhow)
            dire = trace(irhow*b)/rows(irhow)
            idire = totale - dire
        }
        else{
            w = asarray(wx,keys[1])
            totale = sum(irhow*b[1]+irhow*w*b[2])/rows(irhow)
            dire = trace(irhow*b[1]+irhow*w*b[2])/rows(irhow)
            idire = totale - dire
        }
    }
    else{
        NT=0
        totale = 0
        dire = 0
        idire = 0
        for(t=1;t<=T;t++){
            irhow = asarray(iwy,keys[t])
            NT = NT+rows(irhow)
            if (length(b)==1){
                totale = totale + sum(irhow*b)
                dire = dire + trace(irhow*b) 
            }
            else{
                w = asarray(wx,keys[t])
                totale = totale + sum(irhow*b[1]+irhow*w*b[2])
                dire = dire + trace(irhow*b[1]+irhow*w*b[2])
            }
        }
        dire = dire/NT
        totale = totale/NT
        idire = totale - dire

    }
    
    return(totale)

}



// compute [I -rhoW]^(-1) and store in irhow array

void function transinvw(transmorphic matrix wina, 
	                    real scalar T, 
						real scalar rho, 
						transmorphic matrix irhow)
{
	irhow = asarray_create("real")
	keys = asarray_keys(wina)
	if(length(keys)==1){
		w = extrpoi(asarray(wina,keys[1]))
		asarray(irhow,1,matinv(I(rows(w))-rho*w))
	}
	else{
		for(t=1;t<=T;t++){
			w = extrpoi(asarray(wina,keys[t]))
			asarray(irhow,t,matinv(I(rows(w))-rho*w))
		}
	}

}


real colvector rowsd(real matrix x)
{
    var = J(rows(x),1,.)
    for(r=1;r<=rows(x);r++){
        var[r] = sqrt(quadvariance(x[r,.]'))
    }
    return(var)
}

// generate coefficiencts from emprical distribution
real colvector bmc(real matrix V)
{
    return(V*rnormal(rows(V),1,0,1))
}

void function rho_tau(real vector b,
                      real scalar rhoindex,
                      real scalar tauindex,
                      real scalar rhocon,
                      real scalar taucon)
{
    if (rhoindex !=0){
        rhocon = b[rhoindex]
    }
    if (tauindex !=0){
        taucon = b[tauindex]
    }

}

// extract rho,gamma and tau from extimated parameters
void function rho_gamma_tau(real vector b,
                      real scalar rhoindex,
                      real scalar gammaindex,
                      real scalar tauindex,
                      real scalar rhocon,
                      real scalar gammacon,
                      real scalar taucon,
                      real scalar flag)
{
    if(flag!=1){ // flag!=1, no transformation
        if (rhoindex !=0){
            rhocon = b[rhoindex]
        }
        if (tauindex !=0){
            taucon = b[tauindex]
        }
        if (gammaindex !=0){
            gammacon = b[gammaindex]
        }
    }
    else{
        if (rhoindex !=0){
            rhocon = b[rhoindex]
            rymin= st_numscalar("rymin")
            rymax= st_numscalar("rymax")
            rhocon = rymin/(1+exp(rhocon))+rymax*exp(rhocon)/(1+exp(rhocon))
        }
        else{
            rhocon = 0
        }
        if (tauindex !=0){
            taucon = b[tauindex]
            rumin= st_numscalar("rumin")
            rumax= st_numscalar("rumax")
            taucon = rumin/(1+exp(taucon))+rumax*exp(taucon)/(1+exp(taucon))
        }
        else{
            taucon = 0
        }
        if (gammaindex !=0){
            gammacon = b[gammaindex]
            rvmin= st_numscalar("rvmin")
            rvmax= st_numscalar("rvmax")
            gammacon = rvmin/(1+exp(gammacon))+rvmax*exp(gammacon)/(1+exp(gammacon))
        } 
        else{
            gammacon = 0
        }       
    }


}

// loop over specified xvars
real matrix x_mc(string scalar xname,
                     string vector varnames,
                     real colvector bml,
                     real scalar T,
                     transmorphic matrix iwy,
                     transmorphic matrix wxmat,
                     real matrix direxmat,
                     real matrix indirexmat)
    {
        xvar  = tokens(xname)
        totexmat = J(0,1,.)
        direxmat = J(0,1,.)
        indirexmat = J(0,1,.)
        for(i=1;i<=length(xvar);i++){
            ii=extractposx(varnames,xvar[i])
            bi=bml[ii]
            tote=marginx_mc(bi',T,iwy,wxmat,dire=.,indire=.,)
            totexmat=totexmat \ tote
            direxmat=direxmat \ dire
            indirexmat=indirexmat \ indire
        }
        return(totexmat)

    }

    // loop over specified zvars
   real matrix z_mc(string scalar zname,
                     string vector varnames,
                     real colvector bz,
                     real scalar T,
                     real colvector hzust,
                     transmorphic matrix irhow,
                     transmorphic matrix itauw,
                     real matrix direzmat,
                     real matrix indirezmat)
    {
        totezmat = J(0,1,.)
        direzmat = J(0,1,.)
        indirezmat = J(0,1,.)
        zvars = tokens(zname)
        for(i=1;i<=length(zvars);i++){
            ii= select(1..(length(varnames)),varnames:==zvars[i])
            bi=bz[ii]
            tote = marginz_mc(irhow,itauw,bi,hzust,T,dire=.,idire=.,)
            totezmat=totezmat \ tote
            direzmat=direzmat \ dire
            indirezmat=indirezmat \ idire
        }
        return(totezmat)
    }

    real scalar function marginz_mc(transmorphic matrix irhow,
                                    transmorphic matrix itauw,
                                    real scalar bi,
                                    real colvector hzust,
                                    real scalar T,
                                    real scalar dire,
                                    real scalar idire)
    {
        external real vector _pan_tvar
        info = panelsetup(_pan_tvar,1)
        tote = 0
        dire = 0
        idire = 0

        ykeys = asarray_keys(irhow)
        ukeys = asarray_keys(itauw)
        if(length(ykeys)==1){
            iiw1 = asarray(irhow,ykeys[1])
        }
        if(length(ukeys)==1){
            iiw2 = asarray(itauw,ukeys[1])
        }  
        NT = 0  
       
        for(t=1;t<=T;t++){
            ht = panelsubmatrix(hzust,t,info)
            N = length(ht)
            if(length(ykeys)>1){
                iiw1 = asarray(irhow,ykeys[t])
            }
            if(length(ukeys)>1){
                iiw2 = asarray(itauw,ukeys[t])
            }
            tote = tote + sum(iiw1*iiw2*diag(ht)*bi*0.5)
            dire = dire + trace(iiw1*iiw2*diag(ht)*bi*0.5)
            NT = NT + N

        }
        tote = tote/NT
        dire = dire/NT
        idire = tote - dire
        return(tote)
    }



end
