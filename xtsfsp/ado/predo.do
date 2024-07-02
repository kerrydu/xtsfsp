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
local ccomu: list marginvars & uhet   
