#include "Copyright.h"
#include "GAMER.h"

#if ( !defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )




//-----------------------------------------------------------------------------------------
// Function    :  CPU_HydroGravitySolver
// Description :  Use CPU to advance the momentum and energy density by gravitational acceleration
//
// Parameter   :  Flu_Array_New  : Array to store the input and output fluid variables
//                Pot_Array_New  : Array storing the input potential (at the current step)
//                                 --> _New: to be distinguishable from Pot_Array_USG, which is defined at the previous step
//                Corner_Array   : Array storing the physical corner coordinates of each patch
//                Pot_Array_USG  : Array storing the prepared potential          for UNSPLIT_GRAVITY (at the previous step) 
//                Flu_Array_USG  : Array storing the prepared density + momentum for UNSPLIT_GRAVITY (at the previous step)
//                NPatchGroup    : Number of patch groups to be evaluated
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                P5_Gradient    : Use 5-points stencil to evaluate the potential gradient
//                GravityType    : Types of gravity --> self-gravity, external gravity, both
//                ExtAcc_AuxArray: Auxiliary array for adding external acceleration
//                TimeNew        : Physical time at the current  step (for the external gravity solver)
//                TimeOld        : Physical time at the previous step (for the external gravity solver in UNSPLIT_GRAVITY)
//
// Return      :  Flu_Array_New
//-----------------------------------------------------------------------------------------
void CPU_HydroGravitySolver(       real Flu_Array_New[][NCOMP][PS1][PS1][PS1], 
                             const real Pot_Array_New[][GRA_NXT][GRA_NXT][GRA_NXT],
                             const double Corner_Array[][3],
                             const real Pot_Array_USG[][USG_NXT_G][USG_NXT_G][USG_NXT_G], 
                             const real Flu_Array_USG[][GRA_NIN-1][PS1][PS1][PS1],
                             const int NPatchGroup, const real dt, const real dh, const bool P5_Gradient,
                             const OptGravityType_t GravityType, const double ExtAcc_AuxArray[],
                             const double TimeNew, const double TimeOld )
{ 

// check
#  ifdef GAMER_DEBUG
   if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
   if ( TimeNew < 0.0 )
      Aux_Error( ERROR_INFO, "Incorrect TimeNew (%14.7e) !!\n", TimeNew );

#  ifdef UNSPLIT_GRAVITY
   if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
   if ( Pot_Array_USG == NULL || Flu_Array_USG == NULL )
      Aux_Error( ERROR_INFO, "Pot_Array_USG == NULL || Flu_Array_USG == NULL !!\n" );

   if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
   if ( TimeOld >= TimeNew  ||  TimeOld < 0.0 )
      Aux_Error( ERROR_INFO, "Incorrect TimeOld (%14.7e) (TimeNew = %14.7d) !!\n", TimeOld, TimeNew );
#  endif
#  endif


   const int  NPatch    = NPatchGroup*8;
   const real Gra_Const = ( P5_Gradient ) ? -dt/(12.0*dh) : -dt/(2.0*dh);
   const real Const_8   = (real)8.0;

#  ifdef UNSPLIT_GRAVITY
   real AccNew[3], AccOld[3], PxNew, PxOld, PyNew, PyOld, PzNew, PzOld, RhoNew, RhoOld, Eint, _Rho2, Ek;
#  else
   real AccNew[3], PxNew, PyNew, PzNew, RhoNew, Eint, _Rho2;
#  endif
   double x, y, z;


// loop over all patches
#  ifdef UNSPLIT_GRAVITY
#  pragma omp parallel for private( AccNew, AccOld, PxNew, PxOld, PyNew, PyOld, PzNew, PzOld, RhoNew, RhoOld, x, y, z, Eint, _Rho2, Ek ) schedule( runtime )
#  else
#  pragma omp parallel for private( AccNew, PxNew, PyNew, PzNew, RhoNew, x, y, z, Eint, _Rho2 ) schedule( runtime )
#  endif
   for (int P=0; P<NPatch; P++)
   {
#     ifdef UNSPLIT_GRAVITY
      for (int k1=GRA_GHOST_SIZE, k2=USG_GHOST_SIZE, kk=0; k1<GRA_NXT-GRA_GHOST_SIZE; k1++, k2++, kk++)
      for (int j1=GRA_GHOST_SIZE, j2=USG_GHOST_SIZE, jj=0; j1<GRA_NXT-GRA_GHOST_SIZE; j1++, j2++, jj++)
      for (int i1=GRA_GHOST_SIZE, i2=USG_GHOST_SIZE, ii=0; i1<GRA_NXT-GRA_GHOST_SIZE; i1++, i2++, ii++)
#     else
      for (int k1=GRA_GHOST_SIZE, kk=0; k1<GRA_NXT-GRA_GHOST_SIZE; k1++, kk++)
      for (int j1=GRA_GHOST_SIZE, jj=0; j1<GRA_NXT-GRA_GHOST_SIZE; j1++, jj++)
      for (int i1=GRA_GHOST_SIZE, ii=0; i1<GRA_NXT-GRA_GHOST_SIZE; i1++, ii++)
#     endif
      {
         AccNew[0] = (real)0.0;
         AccNew[1] = (real)0.0;
         AccNew[2] = (real)0.0;

#        ifdef UNSPLIT_GRAVITY
         AccOld[0] = (real)0.0;
         AccOld[1] = (real)0.0;
         AccOld[2] = (real)0.0;
#        endif

//       external gravity
         if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
         {
            z = Corner_Array[P][2] + (double)(kk*dh);
            y = Corner_Array[P][1] + (double)(jj*dh);
            x = Corner_Array[P][0] + (double)(ii*dh);

            CPU_ExternalAcc( AccNew, x, y, z, TimeNew, ExtAcc_AuxArray );
            for (int d=0; d<3; d++)    AccNew[d] *= dt;

#           ifdef UNSPLIT_GRAVITY
            CPU_ExternalAcc( AccOld, x, y, z, TimeOld, ExtAcc_AuxArray );
            for (int d=0; d<3; d++)    AccOld[d] *= dt;
#           endif
         }


//       self-gravity
         if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
         {
            if ( P5_Gradient )
            {
               AccNew[0] += Gra_Const * ( -         Pot_Array_New[P][k1  ][j1  ][i1+2] +         Pot_Array_New[P][k1  ][j1  ][i1-2]
                                          + Const_8*Pot_Array_New[P][k1  ][j1  ][i1+1] - Const_8*Pot_Array_New[P][k1  ][j1  ][i1-1] );
               AccNew[1] += Gra_Const * ( -         Pot_Array_New[P][k1  ][j1+2][i1  ] +         Pot_Array_New[P][k1  ][j1-2][i1  ]
                                          + Const_8*Pot_Array_New[P][k1  ][j1+1][i1  ] - Const_8*Pot_Array_New[P][k1  ][j1-1][i1  ] );
               AccNew[2] += Gra_Const * ( -         Pot_Array_New[P][k1+2][j1  ][i1  ] +         Pot_Array_New[P][k1-2][j1  ][i1  ]
                                          + Const_8*Pot_Array_New[P][k1+1][j1  ][i1  ] - Const_8*Pot_Array_New[P][k1-1][j1  ][i1  ] );

#              ifdef UNSPLIT_GRAVITY
               AccOld[0] += Gra_Const * ( -         Pot_Array_USG[P][k2  ][j2  ][i2+2] +         Pot_Array_USG[P][k2  ][j2  ][i2-2]
                                          + Const_8*Pot_Array_USG[P][k2  ][j2  ][i2+1] - Const_8*Pot_Array_USG[P][k2  ][j2  ][i2-1] );
               AccOld[1] += Gra_Const * ( -         Pot_Array_USG[P][k2  ][j2+2][i2  ] +         Pot_Array_USG[P][k2  ][j2-2][i2  ]
                                          + Const_8*Pot_Array_USG[P][k2  ][j2+1][i2  ] - Const_8*Pot_Array_USG[P][k2  ][j2-1][i2  ] );
               AccOld[2] += Gra_Const * ( -         Pot_Array_USG[P][k2+2][j2  ][i2  ] +         Pot_Array_USG[P][k2-2][j2  ][i2  ]
                                          + Const_8*Pot_Array_USG[P][k2+1][j2  ][i2  ] - Const_8*Pot_Array_USG[P][k2-1][j2  ][i2  ] );
#              endif
            }

            else
            {
               AccNew[0] += Gra_Const * ( Pot_Array_New[P][k1  ][j1  ][i1+1] - Pot_Array_New[P][k1  ][j1  ][i1-1] );
               AccNew[1] += Gra_Const * ( Pot_Array_New[P][k1  ][j1+1][i1  ] - Pot_Array_New[P][k1  ][j1-1][i1  ] );
               AccNew[2] += Gra_Const * ( Pot_Array_New[P][k1+1][j1  ][i1  ] - Pot_Array_New[P][k1-1][j1  ][i1  ] );

#              ifdef UNSPLIT_GRAVITY
               AccOld[0] += Gra_Const * ( Pot_Array_USG[P][k2  ][j2  ][i2+1] - Pot_Array_USG[P][k2  ][j2  ][i2-1] );
               AccOld[1] += Gra_Const * ( Pot_Array_USG[P][k2  ][j2+1][i2  ] - Pot_Array_USG[P][k2  ][j2-1][i2  ] );
               AccOld[2] += Gra_Const * ( Pot_Array_USG[P][k2+1][j2  ][i2  ] - Pot_Array_USG[P][k2-1][j2  ][i2  ] );
#              endif
            }
         } // if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )


//       advance fluid
#        ifdef UNSPLIT_GRAVITY
         RhoNew = Flu_Array_New[P][DENS][kk][jj][ii];
         RhoOld = Flu_Array_USG[P][DENS][kk][jj][ii];
         PxNew  = Flu_Array_New[P][MOMX][kk][jj][ii];
         PxOld  = Flu_Array_USG[P][MOMX][kk][jj][ii];
         PyNew  = Flu_Array_New[P][MOMY][kk][jj][ii];
         PyOld  = Flu_Array_USG[P][MOMY][kk][jj][ii];
         PzNew  = Flu_Array_New[P][MOMZ][kk][jj][ii];
         PzOld  = Flu_Array_USG[P][MOMZ][kk][jj][ii];

         _Rho2  = (real)0.5/RhoNew;
         Eint   = Flu_Array_New[P][ENGY][kk][jj][ii] - _Rho2*( SQR(PxNew) + SQR(PyNew) + SQR(PzNew) );

         PxNew += (real)0.5*( RhoOld*AccOld[0] + RhoNew*AccNew[0] );
         PyNew += (real)0.5*( RhoOld*AccOld[1] + RhoNew*AccNew[1] );
         PzNew += (real)0.5*( RhoOld*AccOld[2] + RhoNew*AccNew[2] );

         Flu_Array_New[P][MOMX][kk][jj][ii] = PxNew;
         Flu_Array_New[P][MOMY][kk][jj][ii] = PyNew;
         Flu_Array_New[P][MOMZ][kk][jj][ii] = PzNew;

//       for the unsplitting method, we corret the total energy instead of just kinematic energy (so internal energy may change)
         Flu_Array_New[P][ENGY][kk][jj][ii] += (real)0.5*( PxOld*AccOld[0] + PyOld*AccOld[1] + PzOld*AccOld[2] +
                                                           PxNew*AccNew[0] + PyNew*AccNew[1] + PzNew*AccNew[2] );

//       ensure the positive pressure (restore to the original internal energy if it becomes negative)
         Ek = _Rho2*( SQR(PxNew) + SQR(PyNew) + SQR(PzNew) );
         if ( Flu_Array_New[P][ENGY][kk][jj][ii] - Ek <= (real)0.0 )    Flu_Array_New[P][ENGY][kk][jj][ii] = Eint + Ek;

#        else  // #ifdef UNSPLIT_GRAVITY

         RhoNew = Flu_Array_New[P][DENS][kk][jj][ii];
         PxNew  = Flu_Array_New[P][MOMX][kk][jj][ii];
         PyNew  = Flu_Array_New[P][MOMY][kk][jj][ii];
         PzNew  = Flu_Array_New[P][MOMZ][kk][jj][ii];

         _Rho2  = (real)0.5/RhoNew;
         Eint   = Flu_Array_New[P][ENGY][kk][jj][ii] - _Rho2*( SQR(PxNew) + SQR(PyNew) + SQR(PzNew) );

         PxNew += RhoNew*AccNew[0];
         PyNew += RhoNew*AccNew[1];
         PzNew += RhoNew*AccNew[2];

         Flu_Array_New[P][MOMX][kk][jj][ii] = PxNew;
         Flu_Array_New[P][MOMY][kk][jj][ii] = PyNew;
         Flu_Array_New[P][MOMZ][kk][jj][ii] = PzNew;

//       for the splitting method, we ensure that the internal energy is unchanged
         Flu_Array_New[P][ENGY][kk][jj][ii] = Eint + _Rho2*( SQR(PxNew) + SQR(PyNew) + SQR(PzNew) );
#        endif // #ifdef UNSPLIT_GRAVITY ... else ...

      } // i,j,k
   } // for (int P=0; P<NPatch; P++)

} // FUNCTION : CPU_HydroGravitySolver



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )
