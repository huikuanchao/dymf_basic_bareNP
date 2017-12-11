#include "globals.h"

void bond_stress() {
  int i, j1, j2, m, ind ;
  double mdr2, dr[Dim] ;

  ind = 0 ;

  for (j1=0 ; j1<Dim ; j1++ ) 
    for ( j2=0 ; j2<Dim ; j2++ ){
      Stress_bonds[j1][j2] = 0.0 ;
	for(m=0; m<nthreads; m++)
 	   Stress_bond_t[j1][j2][m] = 0.0 ;
    }

double ratio=double(Nda+Ndb-1)/double(Nda+Ndb);
  // Diblock bonds //
#pragma omp parallel for private(ind,m,mdr2,dr,j1,j2)
  for ( i=0 ; i<nD ; i++ ) {
    for ( m=0 ; m<Nda + Ndb - 1 ; m++ ) {

      ind = i * (Nda + Ndb) + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;
      int tid = omp_get_thread_num() ;
 
      for ( j1=0 ; j1<Dim ; j1++ )
        for ( j2=0 ; j2<Dim ; j2++ )
         Stress_bond_t[j1][j2][tid] += ratio* dr[j1] * dr[j2] ;
      
    }
  
  } // for ( i=0 ; i<nT[k]

 ratio=double(Nha-1)/double(Nha);
  // Homopolymer A bonds //
#pragma omp parallel for private(ind,m,dr,mdr2,j1,j2)
  for ( i=0 ; i<nA ; i++ ) {
    for ( m=0 ; m<Nha - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + i * Nha + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;
	 int tid = omp_get_thread_num() ;
 
      for ( j1=0 ; j1<Dim ; j1++ )
        for ( j2=0 ; j2<Dim ; j2++ )
         Stress_bond_t[j1][j2][tid] += ratio*dr[j1] * dr[j2] ;
      
      
    }
  
  } // for ( i=0 ; i<nT[k]

  // Homopolymer B bonds //
 ratio=double(Nhb-1)/double(Nhb); 
#pragma omp parallel for private(ind,m,mdr2,dr,j1,j2)
  for ( i=0 ; i<nB ; i++ ) {
    for ( m=0 ; m<Nhb - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + nA * Nha + i * Nhb + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;
	 int tid = omp_get_thread_num() ;
 
      for ( j1=0 ; j1<Dim ; j1++ )
        for ( j2=0 ; j2<Dim ; j2++ )
         Stress_bond_t[j1][j2][tid] += ratio*dr[j1] * dr[j2] ;
      
      
    }

  } // for ( i=0 ; i<nT[k]
  
  for ( j1=0 ; j1<Dim ; j1++ )   
    for ( j2=0 ; j2<Dim ; j2++ )
      for ( m=0 ; m<nthreads ; m++ ) 
        Stress_bonds[j1][j2] += 3.0 * Stress_bond_t[j1][j2][m] / V ;



}


void bonds( ) {

  int i, j, k, m, ind ;
  double mdr2, mdr, dr[Dim] ;

  Ubond = 0.0 ;

double ratio  = double(Nda+Ndb-1)/double(Nda+Ndb);
  // Diblock bonds //
#pragma omp parallel for \
  private(ind, j, m, dr, mdr2) \
  reduction(+:Ubond)
  for ( i=0 ; i<nD ; i++ ) {
    for ( m=0 ; m<Nda + Ndb - 1 ; m++ ) {

      ind = i * (Nda + Ndb) + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      Ubond += mdr2 * 1.5 ;

      for ( j=0 ; j<Dim ; j++ ) {
        f[ind][j] -= 3.0 *ratio* dr[j] ;
        f[ind+1][j] += 3.0 * ratio*dr[j] ;
      }

    }

  } // for ( i=0 ; i<nT[k]

ratio = double(Nha-1)/double(Nha);
  // Homopolymer A bonds //
#pragma omp parallel for \
  private(ind, j, m, dr, mdr2) \
  reduction(+:Ubond)
  for ( i=0 ; i<nA ; i++ ) {
    for ( m=0 ; m<Nha - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + i * Nha + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      Ubond += mdr2 * 1.5 ;

      for ( j=0 ; j<Dim ; j++ ) {
        f[ind][j] -= 3.0 * dr[j] ;
        f[ind+1][j] += 3.0 * dr[j] ;
      }

    }

  } // for ( i=0 ; i<nT[k]

  // Homopolymer B bonds //
ratio = double(Nhb-1)/double(Nhb);
#pragma omp parallel for \
  private(ind, j, m, dr, mdr2) \
  reduction(+:Ubond)
  for ( i=0 ; i<nB ; i++ ) {
    for ( m=0 ; m<Nhb - 1 ; m++ ) {

      ind = nD * (Nda + Ndb) + nA * Nha + i * Nhb + m ;

      mdr2 = pbc_mdr2( x[ind] , x[ind+1] , dr ) ;

      Ubond += mdr2 * 1.5 ;

      for ( j=0 ; j<Dim ; j++ ) {
        f[ind][j] -= 3.0 * dr[j] ;
        f[ind+1][j] += 3.0 * dr[j] ;
      }

    }

  } // for ( i=0 ; i<nT[k]


}
