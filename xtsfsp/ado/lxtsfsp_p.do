
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

end