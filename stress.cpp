#include "globals.h"


void calc_stress() {


  int j1, j2 ;
  nonbond_stress() ;
  bond_stress() ;


  for ( j1=0 ; j1<Dim ; j1++ )
    for ( j2=0 ; j2<Dim ; j2++ )
      Ptens[j1][j2] = -(Stress_nb[j1][j2] + Stress_bonds[j1][j2]) + rho0 * Kdelta(j1,j2) ;
      //Ptens[j1][j2] = -Stress_nb[j1][j2]
        //+ double(nA)/V * Kdelta(j1,j2) ;

  Pscalar = 0.0 ;
  for ( j1=0 ; j1<Dim ; j1++ )
    Pscalar += Rg3*Ptens[j1][j1] / double( Dim ) ;

}


void nonbond_stress() {

  int i, j1, j2 ;

  fftw_fwd( rhot , ktmp ) ;
   for(i=0; i<ntypes; i++)
	fftw_fwd( rho[i] , rho_hat[i] ) ;


  for ( j1=0 ; j1<Dim ; j1++ ) {
    for ( j2=0 ; j2<Dim ; j2++ ) {
    //compressibility 
   if(kappa>0){
#pragma omp parallel for  
      for ( i=0 ; i<M ; i++ ){
        ktmp2[i] =  ktmp[i] * vir_func_hat[j1][j2][i] ;
        
      }
      fftw_back( ktmp2 , tmp ) ;

#pragma omp parallel for 
      for ( i=0 ; i<M ; i++ )
        tmp2[i] = rhot[i] * tmp[i]*kappa/2.0/V/rho0 ; 
    
  
     }
// Flory-Huggins
    if(chiAB>0){
#pragma omp parallel for  
      for ( i=0 ; i<M ; i++ ){
        ktmp2[i] = ( rho_hat[1][i]) * vir_func_hat[j1][j2][i] ;
        
      }

      fftw_back( ktmp2 , tmp ) ;
    
#pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        tmp2[i]  += (rho[0][i] )* tmp[i] *chiAB/V/rho0; 
     }

//partilce-particle contribtion 
     if(nP >0 ){      
#pragma omp parallel for  
      for ( i=0 ; i<M ; i++ ){
        ktmp2[i] = ( rho_hat[2][i]) * vir_funcpp_hat[j1][j2][i] ;
      	  
      }

      fftw_back( ktmp2 , tmp ) ;
 #pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        tmp2[i]  += (rho[2][i])* tmp[i] *kappa/V/2.0/rho0; 
     
// partilce-polymer contribution
 #pragma omp parallel for  
      for ( i=0 ; i<M ; i++ ){
        ktmp2[i] = ( rho_hat[2][i]) * vir_funcpg_hat[j1][j2][i] ;
      }

      fftw_back( ktmp2 , tmp ) ;
 #pragma omp parallel for
      for ( i=0 ; i<M ; i++ )
        tmp2[i]  += (rho[0][i]*kappa + (A_partics > 0? chiAB+kappa: kappa)*rho[1][i])* tmp[i] /V/rho0; 
      
     }//nP >0 
      Stress_nb[j1][j2] = -integrate( tmp2 )  ;
    }
  }

}


