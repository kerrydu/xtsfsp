*! Monday, July 1, 2024 at 13:22:19, by Kerry Du
*! total, direct and indirect effects for xtsfsp
capture program drop xtsfsp_margins
program define xtsfsp_margins, rclass

version 16
syntax varlist, [mc(numlist) seed(real 123) NODOTS]

if `"`mc'"'==""{
    xtsfsp_margins_delta `0'
}
else{
    xtsfsp_margins_mc `0'
}
end


capture program drop xtsfsp_margins_delta
program define xtsfsp_margins_delta,rclass
    version 16.0
    syntax varlist 
    local marginvars `varlist'
    ////////////////////////////////////////
    * extract informtaion from estimates
    local cmd  `e(cmd)'
    local cmdline `e(cmdline)'
    local wymat `e(wy)'
    local wxmat `e(wx)'
    local wumat `e(wu)'
    local wxvars `e(wxvars)'
    local T = e(T)
    scalar rymin = e(rymin)
    scalar rymax = e(rymax)
    scalar rumin = e(rumin)
    scalar rumax = e(rumax)
    tempname bml covml
    mat `bml'  = e(b)
    mat `covml' = e(V)
    mata: VCV = st_matrix("`covml'")
    mata: bml = st_matrix("`bml'")
	local varnames: colnames `bml'
    local nparas : word count `varnames'
	mata: varnames=tokens(st_local("varnames"))
    gettoken cmd0 cmdline : cmdline
    cmdlineparse `cmdline' // parse the cmdline of last estimate
    local wy0 `r(wy0)' // original wy information
    local wu0 `r(wu0)' // original wu information
    local wv0 `r(wv0)' // original wv information
    local rhotauindex
    local weightyu
    local taucon = .
    local rhocon = .
    if `"`wu0'"'!= "" {
        local taucon = _b[Wu:_cons]
        local tau_index = `nparas' // para = bx, b_vhet, b_uhet, rho, theta, tau
        local rhotauindex ,`tau_index'
        local weightyu ,`wumat'
    }
    if `"`wy0'"'!= "" {
        local rhocon = _b[Wy:_cons]
        local rho_index = `nparas' // para = bx, b_vhet, b_uhet, rho, theta, tau
        if `"`wu0'"'!= "" {
            local rho_index = `rho_index' -1
        }
        if `"`wv0'"'!= "" {
            local rho_index = `rho_index' -1
        }
        local rhotauindex ,`rho_index' `rhotauindex'
        // if "`wymat'"=="" {
        //     local wymat `wumat' // wy is replaced by wu when wy is not included
        // }
        local weightyu ,`wymat' `weightyu'
    }
    local nx = r(nx)
    local nv = r(nv)
    local nu = r(nu)
    local zvars `r(zvars)'
    local xvars `r(xvars)'
    local nwx: word count `wxvars'
    local nx = `nx' + `nwx'
    *********************
    /////////////////////////////////////////
    if ("`wymat'"=="" & "`wumat'"=="") {
        di as red "No spatial components for the indepvar and inefficiency term are specified in the last estimates" 
        exit
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
        di as red "Check variables specified in frontier and uhet of the last estimates"
		exit 
	}
	local ccomx: list marginvars & xvars
	local ccomu: list marginvars & zvars     

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
   mata: zinuhet = st_data(., "`zvars'")

   if `"`ccomx'"'!=""{
    di "Marginal effects for variables in the frontier: "
    mata: totexmat = J(0,2,.)
    mata: direxmat = J(0,2,.)
    mata: indirexmat = J(0,2,.)

    if ("`wxmat'"==""){
        local wxmat `wymat'
    }

    foreach v in `ccomx'{
        mata: ii=extractposx(varnames[1..`nx'],"`v'")
        mata: bi=bml[ii]
        mata: ii=ii,`rho_index'
        //mata: vi=VCV[ii,ii]
        mata: tote=marginx(`rhocon',bi,VCV[ii,ii],`T',`wymat',`wxmat',dire=.,indire=.,)
        mata: totexmat=totexmat \ tote
        mata: direxmat=direxmat \ dire
        mata: indirexmat=indirexmat \ indire
    }

    mata: totexmat =totexmat, ((totexmat[.,1]):/ totexmat[.,2]), 2*(1:-normal(abs(totexmat[.,1]:/ totexmat[.,2])))
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
    mata: direxmat =direxmat, ((direxmat[.,1]):/ direxmat[.,2]), 2*(1:-normal(abs(direxmat[.,1]:/ direxmat[.,2])))
    mata: st_matrix("direxeff",direxmat)
    mat rownames direxeff = `rnames'
    mat colnames direxeff = "Coeff" "se" "z" "P"			
    dis _n in gr "Direct marginal effect:"
    matlist direxeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
    mata: indirexmat = indirexmat, ((indirexmat[.,1]):/ indirexmat[.,2]), 2*(1:-normal(abs(indirexmat[.,1]:/ indirexmat[.,2])))
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
    di "Marginal effects for variables in the inefficiency function: "
    mata: totezmat = J(0,2,.)
    mata: direzmat = J(0,2,.)
    mata: indirezmat = J(0,2,.)
    local wumat w_ina
    local nxv =`nx' + `nv'
    local nz: word count `zvars'  // nz not include the constant term
    mata: bz=bml[(`nxv'+1)..(`nxv'+`nz')]
    ***->-----<---***
    // three case: 1. rho->(`nxv'+1)..(`nxv'+`nz'), rhoindex ;
    // 2. tau->(`nxv'+1)..(`nxv'+`nz'), tauindex; 
    // 3. rho and tau ->(`nxv'+1)..(`nxv'+`nz'),rhoindex, tauindex
    mata: vii = (`nxv'+1)..(`nxv'+`nz') `rhotauindex' 
    mata: vi = VCV[vii,vii]
    foreach v in `ccomu'{
        mata: ii= select(1..(length(varnames)-`nxv'),varnames[`=`nxv'+1'..length(varnames)]:=="`v'")
        //marginz(`rhocon',`taucon',ii,b,z,vi,`T',dire=.,indire=.,w_ina,`wumat')
        mata: tote = marginz(`rhocon',`taucon',ii,bz,zinuhet,vi,`T',dire=.,indire=. `weightyu')
        mata: totezmat=totezmat \ tote
        mata: direzmat=direzmat \ dire
        mata: indirezmat=indirezmat \ indire
    }
    mata: totezmat =totezmat, ((totezmat[.,1]):/ totezmat[.,2]), 2*(1:-normal(abs(totezmat[.,1]:/ totezmat[.,2])))
    mata: st_matrix("totzeff",totezmat)
    local rnames
    foreach v in `ccomx'{
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
                
    mata: direzmat =direzmat, ((direzmat[.,1]):/ direzmat[.,2]), 2*(1:-normal(abs(direzmat[.,1]:/ direzmat[.,2])))
    mata: st_matrix("direzeff",direzmat)
    mat rownames direzeff = `rnames'
    mat colnames direzeff = "Coeff" "se" "z" "P"			
    dis _n in gr "Direct marginal effect:"
    matlist direzeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
    mata: indirezmat = indirezmat, ((indirezmat[.,1]):/ indirezmat[.,2]), 2*(1:-normal(abs(indirezmat[.,1]:/ indirezmat[.,2])))
    mata: st_matrix("indirezeff",indirezmat)
    mat rownames indirezeff = `rnames'
    mat colnames indirezeff = "Coeff" "se" "z" "P"			
    dis _n in gr "Indirect marginal effect:"
    matlist indirezeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
    return matrix totalzmargins = totzeff
    return matrix directzmargins = direzeff
    return matrix indirectzmargins = indirezeff

}



di  _n "Note: Standard errors are computed using the delta method."
di     "      P value are computed using normal distribution."
restore




end

*************************************************************************************

capture program drop xtsfsp_margins_mc
program define xtsfsp_margins_mc,rclass
version 16.0
    syntax varlist, [mc(real 100) seed(real 123) NODOTS]  
    if "`nodots'"!="" local qui qui
    local marginvars `varlist'
    ////////////////////////////////////////
    * extract informtaion from estimates
    local cmd  `e(cmd)'
    local cmdline `e(cmdline)'
    local wymat `e(wy)'
    local wxmat `e(wx)'
    local wumat `e(wu)'
    local wxvars `e(wxvars)'
    local T = e(T)
    scalar rymin = e(rymin)
    scalar rymax = e(rymax)
    scalar rumin = e(rumin)
    scalar rumax = e(rumax)
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
    local rhotauindex
    local weightyu
    local taucon = .
    local rhocon = .
    local rhotauflag = 0
    local tau_index = 0
    local rho_index = 0
    if `"`wu0'"'!= "" {
        local taucon = _b[Wu:_cons]
        local tau_index = `nparas' // para = bx, b_vhet, b_uhet, rho, theta, tau
        local rhotauindex ,`tau_index'
        local weightyu ,`wumat'
    }
    if `"`wy0'"'!= "" {
        local rhocon = _b[Wy:_cons]
        local rho_index = `nparas' // para = bx, b_vhet, b_uhet, rho, theta, tau
        if `"`wu0'"'!= "" {
            local rho_index = `rho_index' -1
        }
        if `"`wv0'"'!= "" {
            local rho_index = `rho_index' -1
        }
        local rhotauindex ,`rho_index' `rhotauindex'
        // if "`wymat'"=="" {
        //     local wymat `wumat' // wy is replaced by wu when wy is not included
        // }
        local weightyu ,`wymat' `weightyu'
    }
    local nx = r(nx)
    local nv = r(nv)
    local nu = r(nu)
    local zvars `r(zvars)'
    local xvars `r(xvars)'
    local nwx: word count `wxvars'
    local nx = `nx' + `nwx'
    *********************
    /////////////////////////////////////////
    if ("`wymat'"=="" & "`wumat'"=="") {
        di as red "No spatial components for the indepvar and inefficiency term are specified in the last estimates" 
        exit
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
        di as red "Check variables specified in frontier and uhet of the last estimates"
		exit 
	}
	local ccomx: list marginvars & xvars
	local ccomu: list marginvars & zvars     

    confirm  integer number `mc'
    confirm integer number `seed'
    if `mc' < 0{
        di as err "mc() should be a positive integer"
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
   mata: zinuhet = st_data(., "`zvars'")

if `"`ccomx'"'!=""{
    //di "Marginal effects for variables in the frontier: "
    mata: totexmat = J(0,1,.)
    mata: direxmat = J(0,1,.)
    mata: indirexmat = J(0,1,.)

    if ("`wxmat'"==""){
        local wxmat `wymat'
    }

    foreach v in `ccomx'{
        mata: ii=extractposx(varnames[1..`nx'],"`v'")
        mata: bi=bml[ii]
        mata: ii=ii,`rho_index'
        //mata: vi=VCV[ii,ii]
        mata: tote=marginx_mc(`rhocon',bi,`T',`wymat',`wxmat',dire=.,indire=.,)
        mata: totexmat=totexmat \ tote
        mata: direxmat=direxmat \ dire
        mata: indirexmat=indirexmat \ indire
    }

}


if `"`ccomu'"'!=""{
   // di "Marginal effects for variables in the inefficiency function: "
    mata: totezmat = J(0,1,.)
    mata: direzmat = J(0,1,.)
    mata: indirezmat = J(0,1,.)
    local wumat w_ina
    local nxv =`nx' + `nv'
    local nz: word count `zvars'  // nz not include the constant term
    mata: bz=bml[(`nxv'+1)..(`nxv'+`nz')]
 
    foreach v in `ccomu'{
        mata: ii= select(1..(length(varnames)-`nxv'),varnames[`=`nxv'+1'..length(varnames)]:=="`v'")
        //marginz(`rhocon',`taucon',ii,b,z,vi,`T',dire=.,indire=.,w_ina,`wumat')
        mata: tote = marginz_mc(`rhocon',`taucon',ii,bz,zinuhet,`T',dire=.,indire=. `weightyu')
        mata: totezmat=totezmat \ tote
        mata: direzmat=direzmat \ dire
        mata: indirezmat=indirezmat \ indire
    }	//mata: _b_ml = st_matrix("`bml'")
}

** mc to compute standard errors
set seed `seed'
mata: sdxt = J(length(tote),0,.)
mata: sdxd = J(length(tote),0,.)
mata: sdxi = J(length(tote),0,.)
mata: sdzt = J(length(tote),0,.)
mata: sdzd = J(length(tote),0,.)
mata: sdzi = J(length(tote),0,.)
mata: rhocon = .
mata: taucon = .
if `mc'>0 di _n "Monte Carlo simulation: " 
mata: bml0 = bml'
forv b=1/`mc'{
    `qui' displaydot , dot(`b')
    mata: bml =bml0+ bmc(VCV)
    //mata: bml
    mata: rho_tau(bml,`rho_index',`tau_index',rhocon,taucon)
    if `"`ccomx'"'!=""{
        mata: totexb = J(0,1,.)
        mata: direxb = J(0,1,.)
        mata: indirexb = J(0,1,.)
    foreach v in `ccomx'{
        mata: ii=extractposx(varnames[1..`nx'],"`v'")
        mata: bi=bml[ii]
        //mata: vi=VCV[ii,ii]
        mata: tote=marginx_mc(rhocon,bi,`T',`wymat',`wxmat',dire=.,indire=.,)
        mata: totexb=totexb \ tote
        mata: direxb=direxb \ dire
        mata: indirexb=indirexb \ indire
    }
    mata: sdxt = sdxt, totexb
    mata: sdxd = sdxd, direxb
    mata: sdxi = sdxi, indirexb
    }

    if `"`ccomu'"'!=""{
        mata: totezb = J(0,1,.)
        mata: direzb = J(0,1,.)
        mata: indirezb = J(0,1,.)
    foreach v in `ccomu'{
        mata: ii= select(1..(length(varnames)-`nxv'),varnames[`=`nxv'+1'..length(varnames)]:=="`v'")
        //marginz(`rhocon',`taucon',ii,b,z,vi,`T',dire=.,indire=.,w_ina,`wumat')
        mata: tote = marginz_mc(rhocon,taucon,ii,bz,zinuhet,`T',dire=.,indire=. `weightyu')
        mata: totezb=totezb \ tote
        mata: direzb=direzb \ dire
        mata: indirezb=indirezb \ indire
    }
    mata:sdzt = sdzt, totezb
    mata: sdzd = sdzd, direzb
    mata: sdzi = sdzi, indirezb
   }
    
}

 if `"`ccomx'"'!=""{
    mata: sdxt = rowsd(sdxt)
    mata: sdxd = rowsd(sdxd)
    mata: sdxi = rowsd(sdxi)
    mata: totexmat =totexmat,sdxt, (totexmat:/sdxt), 2*(1:-normal(abs(totexmat:/sdxt)))
    mata: direxmat =direxmat, sdxd,(direxmat:/sdxd), 2*(1:-normal(abs(direxmat:/sdxd)))
    mata: indirexmat = indirexmat,sdxi, (indirexmat:/sdxi), 2*(1:-normal(abs(indirexmat:/sdxi)))   
    di "Marginal effects for variables in the frontier: "
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
    di "Marginal effects for variables in the inefficiency function: "
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
	syntax [varname],[NOCONstant]
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


****************************************************************************************

/*!
 * @author Kerry Du
 * @date 2024-07-01
 * @version 1.0
 * @brief This file contains the mata functions for computing the total, direct and indirect effects for xtsfsp
 * @details The total, direct and indirect effects are computed for the frontier and inefficiency function
 * @see xtsfsp_margins_delta, xtsfsp_margins_mc
 */
 
 
cap mata mata drop marginx()
cap mata mata drop marginz()
cap mata mata drop marginx_mc()
cap mata mata drop marginz_mc()
cap mata mata drop transinvw()
cap mata mata drop calceff()
cap mata mata drop calceff_mc()
cap mata mata drop extractposx()
cap mata mata drop rowsd()
cap mata mata drop bmc()
cap mata mata drop rho_tau()
cap mata mata drop IrhoW()
cap mata mata drop IrhoWWIrhoW()
cap mata mata drop IrhoWW()
cap mata mata drop IrhoWWIrhoWW()
mata:
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

// tote=marginx(`rhocon',bi,VCV[ii,ii],`T',w_ina,`wxmat',dire=.,indire=.,)
real vector marginx(real scalar rhocon, 
	                real rowvector b, 
					real matrix V, 
					real scalar T, 
					transmorphic matrix wina, 
					transmorphic matrix wxmat, 
					real vector dire, 
					real vector indire)
{
	IrhoW(wina,T,rhocon,dire0=.,indire0=.) // [I - a*W]^-1
	if(length(b)==1){
		dire = dire0 *b
		indire = indire0*b
		totale = dire + indire
		IrhoWWIrhoW(wina, T, rhocon, dire1=., indire1=.) // [I - a*W]^-1*W*[I - a*W]^-1			
		ddire = dire0 \ dire1
		direse = sqrt(ddire'*V*ddire )
		inddire = indire0 \ indire1
		indirese = sqrt(inddire'*V*inddire )
		totalese = sqrt((ddire+inddire)'*V*(ddire+inddire))		
	}
	else{
		IrhoWW( wina, T, rhocon, dire1=., indire1=.,wxmat)
		dire = dire0 *b[1] + dire1*b[2]
		indire = indire0*b[1] + indire1*b[2]
		totale = dire + indire
		IrhoWWIrhoW(wina, T, rhocon, dire2=., indire2=.)
		IrhoWWIrhoWW(wina, T, rhocon, dire3=., indire3=.,wxmat) // [I - a*W]^-1*W*[I - a*W]^-1*W
		ddire = dire0 \ dire1 \ (dire2*b[1]+dire3*b[2])
		direse = sqrt(ddire'*V*ddire) 
		inddire = indire0 \ indire1 \ (indire2*b[1]+indire3*b[2])
		indirese = sqrt(inddire'*V*inddire)
		totalese =sqrt( (ddire+inddire)'*V*(ddire+inddire))		
	}

	dire = dire, direse
	indire = indire , indirese 
	totale = totale,totalese 
	return(totale)

}


// tote=marginx(`rhocon',bi,VCV[ii,ii],`T',w_ina,`wxmat',dire=.,indire=.,)
real vector marginx_mc(real scalar rhocon, 
                        real rowvector b, 
                        real scalar T, 
                        transmorphic matrix wina, 
                        transmorphic matrix wxmat, 
                        real vector dire, 
                        real vector indire)
{
    IrhoW(wina,T,rhocon,dire0=.,indire0=.) // [I - a*W]^-1
    if(length(b)==1){
        dire = dire0 *b
        indire = indire0*b
        totale = dire + indire	
    }
    else{
        IrhoWW( wina, T, rhocon, dire1=., indire1=.,wxmat)
        dire = dire0 *b[1] + dire1*b[2]
        indire = indire0*b[1] + indire1*b[2]
        totale = dire + indire	
    }

    dire = dire
    indire = indire 
    totale = totale
    return(totale)

}


//marginz(`rhocon',`taucon',ii,b,z,vi,`T',dire=.,indire=.,w_ina,`wumat')


real vector marginz(real scalar rhocon, 
	                real scalar taucon,
					real scalar ii,
					real rowvector b, 
					real matrix z,
					real matrix V, 
					real scalar T,  
					real vector dire, 
					real vector indire,
					transmorphic matrix wymat, |transmorphic matrix wumat)
{

irhow = asarray_create("real")	
if (rhocon ==.){
		asarray(irhow,1,1)
		drhocon = 1
        //wumat = wymat // !!!wy is replaced by wu when wy is not included

}
else{
	rymin = st_numscalar("rymin")
	rymax = st_numscalar("rymax")
	rho = rymin/(1+exp(rhocon))+rymax*exp(rhocon)/(1+exp(rhocon))	
    drhocon =  exp(rhocon)*((rymax-rymin)/(1+exp(rhocon))^2)
	transinvw(wymat,T,rho,irhow)
}
itauw= asarray_create("real")	
if (taucon ==.){
		asarray(itauw,1,1)
		dtaucon =1
}
else{
	if (args()==10){
		rumin = st_numscalar("rymin")
		rumax = st_numscalar("rymax")
		tau = rumin/(1+exp(taucon))+rumax*exp(taucon)/(1+exp(taucon))
		transinvw(wymat,T,tau,itauw)		
	}
	else{
		
		rumin = st_numscalar("rumin")
		rumax = st_numscalar("rumax")
		tau = rumin/(1+exp(taucon))+rumax*exp(taucon)/(1+exp(taucon))
		transinvw(wumat,T,tau,itauw)
	}
		
	dtaucon =  exp(taucon)*((rumax-rumin)/(1+exp(taucon))^2)


}
if (args()==10){
    calceff(irhow,itauw,z,b,ii,T,drhocon,dtaucon,dire=.,indire=.,dvec=.,idvec=.,wymat)
}
else{
    calceff(irhow,itauw,z,b,ii,T,drhocon,dtaucon,dire=.,indire=.,dvec=.,idvec=.,wymat,wumat)
}

dse = sqrt(dvec'*V*dvec)
idse = sqrt(idvec'*V*idvec)
tdse = sqrt((dvec+idvec)'*V*(dvec+idvec))
totale = (dire+indire),tdse
dire = dire, dse
indire = indire , idse 
return(totale)

}




real vector marginz_mc(real scalar rhocon, 
                        real scalar taucon,
                        real scalar ii,
                        real rowvector b, 
                        real matrix z,
                        real scalar T,  
                        real vector dire, 
                        real vector indire,
                        transmorphic matrix wymat, |transmorphic matrix wumat)
{

    irhow = asarray_create("real")	
    if (rhocon ==.){
        asarray(irhow,1,1)
        //wumat = wymat // !!!wy is replaced by wu when wy is not included
    }
    else{
        rymin = st_numscalar("rymin")
        rymax = st_numscalar("rymax")
        rho = rymin/(1+exp(rhocon))+rymax*exp(rhocon)/(1+exp(rhocon))	
        transinvw(wymat,T,rho,irhow)
    }
    itauw= asarray_create("real")	
    if (taucon ==.){
        asarray(itauw,1,1)
    }
    else{
        if (args()==10){
            rumin = st_numscalar("rymin")
            rumax = st_numscalar("rymax")
            tau = rumin/(1+exp(taucon))+rumax*exp(taucon)/(1+exp(taucon))
            transinvw(wymat,T,tau,itauw)		
        }
        else{
            rumin = st_numscalar("rumin")
            rumax = st_numscalar("rumax")
            tau = rumin/(1+exp(taucon))+rumax*exp(taucon)/(1+exp(taucon))
            transinvw(wumat,T,tau,itauw)
        }


    }
    if (args()==10){
        calceff_mc(irhow,itauw,z,b,ii,T,dire=.,indire=.)
    }
    else{
        calceff_mc(irhow,itauw,z,b,ii,T,dire=.,indire=.)
    }

    dire = dire
    indire = indire 
    totale = dire+indire
    return(totale)

}



void function calceff(transmorphic matrix irhow, 
	                  transmorphic matrix itauw, 
					  real matrix z, 
					  real rowvector b, 
					  real scalar ii, 
					  real scalar T, 
					  real scalar drhocon, 
					  real scalar dtaucon,
					  real vector dire, 
					  real vector indire,
					  real vector dvec,
					  real vector idvec,
					  transmorphic matrix wyina,| transmorphic matrix wumat)
{
	external real colvector _pan_tvar
	info = panelsetup(_pan_tvar,1)
	//nt = panelstats(info)[1]
	wykeys = asarray_keys(irhow)
	wukeys = asarray_keys(itauw)
	NN=0
	dire = 0
	indire = 0
	drho = 0
	idrho = 0
	dtau =0
	idtau = 0
	ddelta = 0
	iddelta = 0
	for(t=1;t<=T;t++){
		zt = panelsubmatrix(z,t,info)
		if(length(wykeys)==1){
			iiw1 = asarray(irhow,wykeys[1])
			wy = asarray(wyina,wykeys[1])
		}
		else{
			iiw1 = asarray(irhow,wykeys[t])
			wy = asarray(wyina,wykeys[t])
		}

		if(length(wukeys)==1){
			iiw2 =asarray(itauw,wukeys[1])
			if (args()==13){
				wu = wy
			}
			else{
				wu = asarray(wumat,wukeys[1])
			}
		}
		else{
			iiw2 = asarray(itauw,wukeys[t])
			if (args()==13){
				wu = wy
			}
			else{
				wu = asarray(wumat,wukeys[t])
			}
		}
		iiw = iiw1*iiw2
		ht = exp(zt*b'*0.5)
		N = rows(zt)
		NN = NN + N
		dire = dire + trace(iiw*diag(ht))*0.5*b[ii]
		indire = indire +(sum(iiw*diag(ht)) - trace(iiw*diag(ht)))*0.5*b[ii]
		if (rows(iiw1)!=1){
			d0 = iiw1*wy*iiw1*iiw2*diag(ht)*0.5*b[ii]*drhocon
			drho = drho + trace(d0)
			idrho = idrho + sum(d0) - trace(d0)
		}
		if(rows(iiw2)!=1){
			d0 = iiw1*iiw2*wu*iiw2*diag(ht)*0.5*b[ii]*dtaucon
			dtau = dtau + trace(d0)
			idtau = idtau + sum(d0) - trace(d0)
		}
		ddelta = J(cols(z),1,0)
		iddelta = J(cols(z),1,0)
		for(j=1;j<=cols(z);j++){
			d0 = iiw1*iiw2*diag(ht:*zt[.,j]*0.5)*0.5*b[ii]
			if(j==ii){
				d0 = d0 + iiw1*iiw2*diag(ht)*0.5
			}
			ddelta[j] =ddelta[j]+ trace(d0)
			iddelta[j] =iddelta[j]+ sum(d0) - trace(d0)
		}
	}
    ddelta = ddelta/NN
    iddelta = iddelta/NN
    drho = drho/NN
    idrho = idrho/NN
    dtau = dtau/NN
    idtau = idtau/NN
	dire = dire/NN 
	indire = indire /NN
	if (rows(iiw1)!=1){
		dvec = ddelta\ drho
		idvec = iddelta\ idrho
	}
	if (rows(iiw2)!=1){
		dvec = dvec \  dtau
		idvec = idvec\  idtau
	}

}


void function calceff_mc(transmorphic matrix irhow, 
                        transmorphic matrix itauw, 
                        real matrix z, 
                        real rowvector b, 
                        real scalar ii, 
                        real scalar T, 
                        real vector dire, 
                        real vector indire)
{
    external real colvector _pan_tvar
    info = panelsetup(_pan_tvar,1)
    //nt = panelstats(info)[1]
    wykeys = asarray_keys(irhow)
    wukeys = asarray_keys(itauw)
    NN=0
    dire = 0
    indire = 0
for(t=1;t<=T;t++){
    zt = panelsubmatrix(z,t,info)
    if(length(wykeys)==1){
        iiw1 = asarray(irhow,wykeys[1])
    }
    else{
        iiw1 = asarray(irhow,wykeys[t])
    }

    if(length(wukeys)==1){
        iiw2 =asarray(itauw,wukeys[1])
    }
    else{
        iiw2 = asarray(itauw,wukeys[t])
    }
    iiw = iiw1*iiw2
    ht = exp(zt*b'*0.5)
    N = rows(zt)
    NN = NN + N
    dire = dire + trace(iiw*diag(ht))*0.5*b[ii]
    indire = indire +(sum(iiw*diag(ht)) - trace(iiw*diag(ht)))*0.5*b[ii]

}
dire = dire/NN 
indire = indire /NN

}


void function transinvw(transmorphic matrix wina, 
	                    real scalar T, 
						real scalar rho, 
						transmorphic matrix irhow)
{
	//irhow = asarray_create("real")
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

//old codes

//////////////////
// rho_con=>rho 
// 对rho求导=》rho_con 
// 求导向量为(beta;theta;rho)

// 对beta求导
void function IrhoW( transmorphic matrix wina,real scalar T, ///
    real scalar rhocon,real scalar dire, real scalar indire)
{
rmin = st_numscalar("rmin")
rmax = st_numscalar("rmax")	
rho = rmin/(1+exp(rhocon))+rmax*exp(rhocon)/(1+exp(rhocon))

wkeys = sort(asarray_keys(wina),1)

if(length(wkeys)==1){

w=extrpoi(asarray(wina,wkeys[1]))
N=rows(w)
wirho = matinv(I(N)-rho*w)
dire = trace(wirho)/N
indire = (sum(wirho) - trace(wirho))/N
}
else{
NN=0
//NNN=0
dire=0
indire=0
for(j=1;j<=T;j++){

w=extrpoi(asarray(wina,wkeys[j]))
N=rows(w)
NN = NN +N
wirho = matinv(I(N)-rho*w)
dire =dire+ trace(wirho)
indire =indire+ (sum(wirho) - trace(wirho))
//NNN = NNN + N*(N-1)

}
dire = dire/NN 
indire = indire /NN
}

}

//对theta求导
void function IrhoWW( transmorphic matrix wina,real scalar T, ///
    real scalar rhocon,real scalar dire, ///
    real scalar indire, |  transmorphic matrix wx)
{

rmin = st_numscalar("rmin")
rmax = st_numscalar("rmax")	
rho = rmin/(1+exp(rhocon))+rmax*exp(rhocon)/(1+exp(rhocon))	

if(args()==5){
wkeys = sort(asarray_keys(wina),1)
if(length(wkeys)==1){

w=extrpoi(asarray(wina,wkeys[1]))
N=rows(w)
wirho = matinv(I(N)-rho*w)*w
dire = trace(wirho)/N
indire = (sum(wirho) - trace(wirho))/N
}
else{
NN=0
//NNN=0
dire=0
indire=0
for(j=1;j<=T;j++){

w=extrpoi(asarray(wina,wkeys[j]))
N=rows(w)
NN = NN +N
wirho = matinv(I(N)-rho*w)*w
dire =dire+ trace(wirho)
indire =indire+ (sum(wirho) - trace(wirho))
// NNN = NNN + N*(N-1)

}
dire = dire/NN 
indire = indire /NN
}		

}
else{

NN=0
//NNN=0
dire=0
indire=0
wykeys = sort(asarray_keys(wina),1)
wxkeys = sort(asarray_keys(wx),1)
for(j=1;j<=T;j++){
if(length(wykeys)==1) w=extrpoi(asarray(wina,wykeys[1]))
else w=extrpoi(asarray(wina,wykeys[j]))
if (length(wxkeys)==1) w2=extrpoi(asarray(wx,wxkeys[1]))
else w2=extrpoi(asarray(wx,wxkeys[j]))
N=rows(w)
NN = NN +N
wirho = matinv(I(N)-rho*w)*w2
dire =dire+ trace(wirho)
indire =indire+ (sum(wirho) - trace(wirho))
//NNN = NNN + N*(N-1)
}
dire = dire/NN 
indire = indire /NN		

}



}



//对rho求导1
void function IrhoWWIrhoW( transmorphic matrix wina,real scalar T, ///
    real scalar rhocon,real scalar dire, real scalar indire)
{
rmin = st_numscalar("rmin")
rmax = st_numscalar("rmax")	
rho = rmin/(1+exp(rhocon))+rmax*exp(rhocon)/(1+exp(rhocon))
drho =  exp(rhocon)*((rmax-rmin)/(1+exp(rhocon))^2)
wkeys = sort(asarray_keys(wina),1)

if(length(wkeys)==1){

w=extrpoi(asarray(wina,wkeys[1]))
N=rows(w)
wirho = matinv(I(N)-rho*w)*w*matinv(I(N)-rho*w)*drho
dire = trace(wirho)/N
indire = (sum(wirho) - trace(wirho))/N
}
else{
NN=0
//NNN=0
dire=0
indire=0
for(j=1;j<=T;j++){

w=extrpoi(asarray(wina,wkeys[j]))
N=rows(w)
NN = NN +N
wirho = matinv(I(N)-rho*w)*w*matinv(I(N)-rho*w)*drho
dire =dire+ trace(wirho)
indire =indire+ (sum(wirho) - trace(wirho))
//NNN = NNN + N*(N-1)

}
dire = dire/NN 
indire = indire /NN
}

}


//对rho求导2
void function IrhoWWIrhoWW( transmorphic matrix wina,real scalar T, ///
    real scalar rhocon,real scalar dire, real scalar indire, | transmorphic matrix wx)
{

rmin = st_numscalar("rmin")
rmax = st_numscalar("rmax")	
rho = rmin/(1+exp(rhocon))+rmax*exp(rhocon)/(1+exp(rhocon))	
drho =  exp(rhocon)*((rmax-rmin)/(1+exp(rhocon))^2)
if (args()==5){

wkeys = sort(asarray_keys(wina),1)

if(length(wkeys)==1){

w=extrpoi(asarray(wina,wkeys[1]))
N=rows(w)
wirho = matinv(I(N)-rho*w)*w*matinv(I(N)-rho*w)*w*drho
dire = trace(wirho)/N
indire = (sum(wirho) - trace(wirho))/N
}
else{
NN=0
//NNN=0
dire=0
indire=0
for(j=1;j<=T;j++){

w=extrpoi(asarray(wina,wkeys[j]))
N=rows(w)
NN = NN +N
wirho = matinv(I(N)-rho*w)*w*matinv(I(N)-rho*w)*w*drho
dire =dire+ trace(wirho)
indire =indire+ (sum(wirho) - trace(wirho))
//NNN = NNN + N*(N-1)

}
dire = dire/NN 
indire = indire /NN
}		

}
else{

NN=0
//NNN=0
dire=0
indire=0
wykeys = sort(asarray_keys(wina),1)
wxkeys = sort(asarray_keys(wx),1)
for(j=1;j<=T;j++){
if(length(wykeys)==1) w=extrpoi(asarray(wina,wykeys[1]))
else w=extrpoi(asarray(wina,wykeys[j]))
if (length(wxkeys)==1) w2=extrpoi(asarray(wx,wxkeys[1]))
else w2=extrpoi(asarray(wx,wxkeys[j]))
N=rows(w)
NN = NN +N
wirho = matinv(I(N)-rho*w)*w*matinv(I(N)-rho*w)*w2*drho
dire =dire+ trace(wirho)
indire =indire+ (sum(wirho) - trace(wirho))
//NNN = NNN + N*(N-1)
}
dire = dire/NN 
indire = indire /NN

}



}


end