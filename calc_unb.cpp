#include "globals.h"

void calc_Unb() {

  int i ;
  // Take Gaussian potential to k-space
  fftw_fwd( uG , ktmp2 ) ;
  
  
  // Polymer chi A-B contribution //
  fftw_fwd( rho[0] , ktmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    ktmp[i] *= ktmp2[i] ;

  fftw_back( ktmp , tmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] *= rho[1][i] ;

  U_chi_gg = integrate( tmp ) * chiAB / rho0 ;


  // Polymer kappa contribution //
#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp2[i] = rhot[i] ;//- rho0 ;

  fftw_fwd( tmp2 , ktmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    ktmp[i] *= ktmp2[i] ;

  fftw_back( ktmp , tmp ) ;

#pragma omp parallel for
  for ( i=0 ; i<M ; i++ )
    tmp[i] *= tmp2[i] ;

  U_kappa_gg = integrate( tmp ) * kappa / 2.0 / rho0 ;


// Nparticle-Nparticle contribution // 

  if( nP > 0 ){
        fftw_fwd(uP, ktmp2);
	fftw_fwd( rho[2] , ktmp ) ;
 #pragma omp parallel for
 	for ( i=0 ; i<M ; i++ )       
	   ktmp2[i] *= ktmp[i] ;
        fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for
	for ( i=0 ; i<M ; i++ )
	     tmp[i] *= rho[2][i];
	
	U_kappa_pp = integrate( tmp ) * kappa / 2.0 /rho0 ;
        //A like particle ?//
  
  	fftw_fwd(uPG, ktmp2);
#pragma omp parallel for
		for ( i=0 ; i<M ; i++ )
  		   ktmp2[i] *= ktmp[i] ;
        	fftw_back( ktmp2 , tmp ) ;
#pragma omp parallel for
		for ( i=0 ; i<M ; i++ )
	     		tmp[i] *= ( (A_partics >0 ? chiAB: 0)+kappa/2.0)*rho[1][i]/rho0 + kappa*rho[0][i]/rho0;
	
	U_kappa_pg  = integrate( tmp )  ;
 	    
  
  }



}

