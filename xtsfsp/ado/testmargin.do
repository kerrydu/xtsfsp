
mat b =e(b)
mat list b
mat V = e(V)
mat list V
local varnames: colnames b 
di "`varnames'"

local marginvars x 



//do predo.do 
do lxtsfsp_p.do
local marginvars z
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
di "`cmdline'"
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
    di "----`weightyu'----"
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
    di "----`weightyu'----"
}

di "`weightyu'"



local nx = r(nx)
local nv = r(nv)
local nu = r(nu)
local zvars `r(zvars)'
di "`zvars'"
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
   di "...`zvars'..."
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
        mata: tote=marginx(`rhocon',bi,VCV[ii,ii],`T',`wymat',`wxmat',dire=.,indire=.)
        mata: totexmat=totexmat \ tote
        mata: direxmat=direxmat \ dire
        mata: indirexmat=indirexmat \ indire
    }

    mata: totexmat =totexmat, ((totexmat[.,1]):/ totexmat[.,2]), 2*(1-normal(abs(totexmat[.,1]:/ totexmat[.,2])))
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
    mata: direxmat =direxmat, ((direxmat[.,1]):/ direxmat[.,2]), 2*(1-normal(abs(direxmat[.,1]:/ direxmat[.,2])))
    mata: st_matrix("direxeff",direxmat)
    mat rownames direxeff = `rnames'
    mat colnames direxeff = "Coeff" "se" "z" "P"			
    dis _n in gr "Direct marginal effect:"
    matlist direxeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
    mata: indirexmat = indirexmat, ((indirexmat[.,1]):/ indirexmat[.,2]), 2*(1-normal(abs(indirexmat[.,1]:/ indirexmat[.,2])))
    mata: st_matrix("indirexeff",indirexmat)
    mat rownames indirexeff = `rnames'
    mat colnames indirexeff = "Coeff" "se" "z" "P"		
    dis _n in gr "Indirect marginal effect:"
    matlist indirexeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")

}



if `"`ccomu'"'!=""{
    di "Marginal effects for variables in the inefficiency function: "
    mata: totezmat = J(0,2,.)
    mata: direzmat = J(0,2,.)
    mata: indirezmat = J(0,2,.)
    local wumat w_ina
    local nxv =`nx' + `nv'
    local nz: word count `zvars'  // nz not include the constant term
    di "(`nxv'+1)..(`nxv'+`nz')"
    mata: bml
    mata: bz=bml[(`nxv'+1)..(`nxv'+`nz')]
    ***->-----<---***
    // three case: 1. rho->(`nxv'+1)..(`nxv'+`nz'), rhoindex ;
    // 2. tau->(`nxv'+1)..(`nxv'+`nz'), tauindex; 
    // 3. rho and tau ->(`nxv'+1)..(`nxv'+`nz'),rhoindex, tauindex
    mata: vii = (`nxv'+1)..(`nxv'+`nz') `rhotauindex' 
    mata: vii
    mata: vi = VCV[vii,vii ]
    mata: vi
    foreach v in `ccomu'{
        mata: ii= select(1..(length(varnames)-`nxv'),varnames[`=`nxv'+1'..length(varnames)]:=="`v'")
        //marginz(`rhocon',`taucon',ii,b,z,vi,`T',dire=.,indire=.,w_ina,`wumat')
        mata: vi
        mata: tote = marginz(`rhocon',`taucon',ii,bz,zinuhet,vi,`T',dire=.,indire=. `weightyu')
        mata: rows(tote),cols(tote)
        mata: totezmat=totezmat \ tote
        mata: direzmat=direzmat \ dire
        mata: indirezmat=indirezmat \ indire
    }
    mata: totezmat 
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
                
    mata: direzmat =direzmat, ((direzmat[.,1]):/ direzmat[.,2]), 2*(1-normal(abs(direzmat[.,1]:/ direzmat[.,2])))
    mata: st_matrix("direzeff",direzmat)
    mat rownames direzeff = `rnames'
    mat colnames direzeff = "Coeff" "se" "z" "P"			
    dis _n in gr "Direct marginal effect:"
    matlist direzeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")
    mata: indirezmat = indirezmat, ((indirezmat[.,1]):/ indirezmat[.,2]), 2*(1-normal(abs(indirezmat[.,1]:/ indirezmat[.,2])))
    mata: st_matrix("indirezeff",indirezmat)
    mat rownames indirezeff = `rnames'
    mat colnames indirezeff = "Coeff" "se" "z" "P"			
    dis _n in gr "Indirect marginal effect:"
    matlist indirezeff, cspec(`cf') rspec(`rf') noblank rowtitle("Variable")

}

