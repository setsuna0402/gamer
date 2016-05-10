#include "Copyright.h"
#include "Macro.h"
#include "CUPOT.h"

#if ( defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )



#include "../../SelfGravity/GPU_Gravity/CUPOT_ExternalAcc.cu"
#define GRA_NTHREAD  ( PATCH_SIZE*PATCH_SIZE*GRA_BLOCK_SIZE_Z )

// variables reside in constant memory
__constant__ double ExtAcc_AuxArray_d[EXT_ACC_NAUX_MAX];




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_HydroGravitySolver_SetConstMem
// Description :  Set the constant memory used by CUPOT_HydroGravitySolver
//
// Note        :  Adopt the suggested approach for CUDA version >= 5.0
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
int CUPOT_HydroGravitySolver_SetConstMem( double ExtAcc_AuxArray_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( ExtAcc_AuxArray_d, ExtAcc_AuxArray_h, EXT_ACC_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUPOT_HydroGravitySolver_SetConstMem



//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_HydroGravitySolver
// Description :  GPU gravity solver which advances the momentum and energy density of a group of patches 
//                by gravitational acceleration (including the external gravity)
//
// Note        :  Prefix "g" for pointers pointing to the "Global" memory space
//                Prefix "s" for pointers pointing to the "Shared" memory space
//
// Parameter   :  g_Flu_Array_New   : Global memory array to store the input and output fluid variables
//                g_Pot_Array_New   : Global memory array storing the input potential for evaluating the 
//                                    gravitational acceleration
//                g_Corner_Array    : Global memory array storing the physical corner coordinates of each patch
//                g_Pot_Array_USG   : Array storing the prepared potential          for UNSPLIT_GRAVITY (at the previous step) 
//                g_Flu_Array_USG   : Array storing the prepared density + momentum for UNSPLIT_GRAVITY (at the previous step)
//                Gra_Const         : 3-P stencil : -dt / ( 2*dh) 
//                                    5-P stencil : -dt / (12*dh)
//                P5_Gradient       : Use 5-points stecil to evaluate the potential gradient
//                GravityType       : Types of gravity --> self-gravity, external gravity, both
//                TimeNew           : Physical time at the current  step (for the external gravity solver)
//                TimeOld           : Physical time at the previous step (for the external gravity solver in UNSPLIT_GRAVITY)
//                dt                : Evolution time-step (for the external gravity solver)
//                dh                : Grid size (for the external gravity solver)
//---------------------------------------------------------------------------------------------------
__global__ void CUPOT_HydroGravitySolver(       real g_Flu_Array_New[][NCOMP][ PATCH_SIZE*PATCH_SIZE*PATCH_SIZE ],
                                          const real g_Pot_Array_New[][ GRA_NXT*GRA_NXT*GRA_NXT ],
                                          const double g_Corner_Array[][3],
                                          const real g_Pot_Array_USG[][ USG_NXT_G*USG_NXT_G*USG_NXT_G ],
                                          const real g_Flu_Array_USG[][GRA_NIN-1][ PS1*PS1*PS1 ],
                                          const real Gra_Const, const bool P5_Gradient, const OptGravityType_t GravityType, 
                                          const double TimeNew, const double TimeOld, const real dt, const real dh )
{

   const uint bx     = blockIdx.x;
   const uint tx     = threadIdx.x; 
   const uint ty     = threadIdx.y; 
   const uint tz     = threadIdx.z; 
   const uint ID     = __umul24( tz, PS1*PS1 ) + __umul24( ty, PS1 ) + tx;
   const uint NSlice = GRA_BLOCK_SIZE_Z;

   uint g_idx        = ID; 
   uint s_idx_new    =   __umul24( GRA_GHOST_SIZE+tz, GRA_NXT*GRA_NXT ) 
                       + __umul24( GRA_GHOST_SIZE+ty, GRA_NXT   ) + (GRA_GHOST_SIZE+tx);

   uint   ip1_new, jp1_new, kp1_new, im1_new, jm1_new, km1_new, t;
   uint   ip2_new, jp2_new, kp2_new, im2_new, jm2_new, km2_new;
   real   Acc_new[3], px_new, py_new, pz_new, rho_new, Eint, _rho2;
   double x, y, z;

   __shared__ real s_Pot_new[GRA_NXT*GRA_NXT*GRA_NXT];

#  ifdef UNSPLIT_GRAVITY
   uint s_idx_old    =   __umul24( USG_GHOST_SIZE+tz, USG_NXT_G*USG_NXT_G ) 
                       + __umul24( USG_GHOST_SIZE+ty, USG_NXT_G ) + (USG_GHOST_SIZE+tx);

   uint ip1_old, jp1_old, kp1_old, im1_old, jm1_old, km1_old;
   uint ip2_old, jp2_old, kp2_old, im2_old, jm2_old, km2_old;
   real Acc_old[3], px_old, py_old, pz_old, rho_old, Ek;

   __shared__ real s_Pot_old[USG_NXT_G*USG_NXT_G*USG_NXT_G];
#  endif


// set the physical coordinates of each cell for the external gravity solver
   if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
   {
      x = g_Corner_Array[bx][0] + (double)(tx*dh);
      y = g_Corner_Array[bx][1] + (double)(ty*dh);
      z = g_Corner_Array[bx][2] + (double)(tz*dh);
   }


// load the potential from the global memory to the shared memory 
   if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
   {
      t = ID;
      do {  s_Pot_new[t] = g_Pot_Array_New[bx][t];   t += GRA_NTHREAD; }  while ( t < CUBE(GRA_NXT) );

      t = ID;
#     ifdef UNSPLIT_GRAVITY
      do {  s_Pot_old[t] = g_Pot_Array_USG[bx][t];   t += GRA_NTHREAD; }  while ( t < CUBE(USG_NXT_G) );
#     endif
   }

   __syncthreads();

   
   for (uint Slice=tz; Slice<PS1; Slice+=NSlice)
   {
      ip1_new = s_idx_new + 1;
      jp1_new = s_idx_new + GRA_NXT;
      kp1_new = s_idx_new + GRA_NXT*GRA_NXT;
      im1_new = s_idx_new - 1;
      jm1_new = s_idx_new - GRA_NXT;
      km1_new = s_idx_new - GRA_NXT*GRA_NXT;

#     ifdef UNSPLIT_GRAVITY
      ip1_old = s_idx_old + 1;
      jp1_old = s_idx_old + USG_NXT_G;
      kp1_old = s_idx_old + USG_NXT_G*USG_NXT_G;
      im1_old = s_idx_old - 1;
      jm1_old = s_idx_old - USG_NXT_G;
      km1_old = s_idx_old - USG_NXT_G*USG_NXT_G;
#     endif

      if ( P5_Gradient )
      {
         ip2_new = s_idx_new + 2;
         jp2_new = s_idx_new + 2*GRA_NXT;
         kp2_new = s_idx_new + 2*GRA_NXT*GRA_NXT;
         im2_new = s_idx_new - 2;
         jm2_new = s_idx_new - 2*GRA_NXT;
         km2_new = s_idx_new - 2*GRA_NXT*GRA_NXT;

#        ifdef UNSPLIT_GRAVITY
         ip2_old = s_idx_old + 2;
         jp2_old = s_idx_old + 2*USG_NXT_G;
         kp2_old = s_idx_old + 2*USG_NXT_G*USG_NXT_G;
         im2_old = s_idx_old - 2;
         jm2_old = s_idx_old - 2*USG_NXT_G;
         km2_old = s_idx_old - 2*USG_NXT_G*USG_NXT_G;
#        endif
      } // if ( P5_Gradient )

      
//    1. evalute the gravitational acceleration
      Acc_new[0] = (real)0.0;
      Acc_new[1] = (real)0.0;
      Acc_new[2] = (real)0.0;

#     ifdef UNSPLIT_GRAVITY
      Acc_old[0] = (real)0.0;
      Acc_old[1] = (real)0.0;
      Acc_old[2] = (real)0.0;
#     endif

//    1.1 external gravity
      if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
      {
         CUPOT_ExternalAcc( Acc_new, x, y, z, TimeNew, ExtAcc_AuxArray_d );
         for (int d=0; d<3; d++)    Acc_new[d] *= dt;

#        ifdef UNSPLIT_GRAVITY
         CUPOT_ExternalAcc( Acc_old, x, y, z, TimeOld, ExtAcc_AuxArray_d );
         for (int d=0; d<3; d++)    Acc_old[d] *= dt;
#        endif
      }

//    1.2 self-gravity
      if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
      {
         if ( P5_Gradient )   // 5-point gradient
         {
            Acc_new[0] += Gra_Const*( - s_Pot_new[ip2_new] + (real)8.0*s_Pot_new[ip1_new] - (real)8.0*s_Pot_new[im1_new] + s_Pot_new[im2_new] );
            Acc_new[1] += Gra_Const*( - s_Pot_new[jp2_new] + (real)8.0*s_Pot_new[jp1_new] - (real)8.0*s_Pot_new[jm1_new] + s_Pot_new[jm2_new] );
            Acc_new[2] += Gra_Const*( - s_Pot_new[kp2_new] + (real)8.0*s_Pot_new[kp1_new] - (real)8.0*s_Pot_new[km1_new] + s_Pot_new[km2_new] );

#           ifdef UNSPLIT_GRAVITY
            Acc_old[0] += Gra_Const*( - s_Pot_old[ip2_old] + (real)8.0*s_Pot_old[ip1_old] - (real)8.0*s_Pot_old[im1_old] + s_Pot_old[im2_old] );
            Acc_old[1] += Gra_Const*( - s_Pot_old[jp2_old] + (real)8.0*s_Pot_old[jp1_old] - (real)8.0*s_Pot_old[jm1_old] + s_Pot_old[jm2_old] );
            Acc_old[2] += Gra_Const*( - s_Pot_old[kp2_old] + (real)8.0*s_Pot_old[kp1_old] - (real)8.0*s_Pot_old[km1_old] + s_Pot_old[km2_old] );
#           endif
         }

         else                 // 3-point gradient
         {
            Acc_new[0] += Gra_Const*( s_Pot_new[ip1_new] - s_Pot_new[im1_new] );
            Acc_new[1] += Gra_Const*( s_Pot_new[jp1_new] - s_Pot_new[jm1_new] );
            Acc_new[2] += Gra_Const*( s_Pot_new[kp1_new] - s_Pot_new[km1_new] );

#           ifdef UNSPLIT_GRAVITY
            Acc_old[0] += Gra_Const*( s_Pot_old[ip1_old] - s_Pot_old[im1_old] );
            Acc_old[1] += Gra_Const*( s_Pot_old[jp1_old] - s_Pot_old[jm1_old] );
            Acc_old[2] += Gra_Const*( s_Pot_old[kp1_old] - s_Pot_old[km1_old] );
#           endif
         }
      } // if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )


//    2. advance the fluid
#     ifdef UNSPLIT_GRAVITY
      rho_new = g_Flu_Array_New[bx][DENS][g_idx];
      rho_old = g_Flu_Array_USG[bx][DENS][g_idx];
      px_new  = g_Flu_Array_New[bx][MOMX][g_idx];
      px_old  = g_Flu_Array_USG[bx][MOMX][g_idx];
      py_new  = g_Flu_Array_New[bx][MOMY][g_idx];
      py_old  = g_Flu_Array_USG[bx][MOMY][g_idx];
      pz_new  = g_Flu_Array_New[bx][MOMZ][g_idx];
      pz_old  = g_Flu_Array_USG[bx][MOMZ][g_idx];

      _rho2   = (real)0.5/rho_new;
      Eint    = g_Flu_Array_New[bx][ENGY][g_idx] - _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );

      px_new += (real)0.5*( rho_old*Acc_old[0] + rho_new*Acc_new[0] );
      py_new += (real)0.5*( rho_old*Acc_old[1] + rho_new*Acc_new[1] );
      pz_new += (real)0.5*( rho_old*Acc_old[2] + rho_new*Acc_new[2] );

      g_Flu_Array_New[bx][MOMX][g_idx] = px_new;
      g_Flu_Array_New[bx][MOMY][g_idx] = py_new;
      g_Flu_Array_New[bx][MOMZ][g_idx] = pz_new;

//    for the unsplitting method, we corret the total energy instead of just kinematic energy (so internal energy may change)
      g_Flu_Array_New[bx][ENGY][g_idx] += (real)0.5*( px_old*Acc_old[0] + py_old*Acc_old[1] + pz_old*Acc_old[2] +
                                                      px_new*Acc_new[0] + py_new*Acc_new[1] + pz_new*Acc_new[2] );

//    ensure the positive pressure (restore to the original internal energy if it becomes negative)
      Ek = _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );
      if ( g_Flu_Array_New[bx][ENGY][g_idx] - Ek <= (real)0.0 )   g_Flu_Array_New[bx][ENGY][g_idx] = Eint + Ek;

#     else // #ifdef UNSPLIT_GRAVITY

      rho_new = g_Flu_Array_New[bx][DENS][g_idx];
      px_new  = g_Flu_Array_New[bx][MOMX][g_idx];
      py_new  = g_Flu_Array_New[bx][MOMY][g_idx];
      pz_new  = g_Flu_Array_New[bx][MOMZ][g_idx];

      _rho2   = (real)0.5/rho_new;
      Eint    = g_Flu_Array_New[bx][ENGY][g_idx] - _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );
      
      px_new += rho_new*Acc_new[0];
      py_new += rho_new*Acc_new[1];
      pz_new += rho_new*Acc_new[2];
      
      g_Flu_Array_New[bx][MOMX][g_idx] = px_new;
      g_Flu_Array_New[bx][MOMY][g_idx] = py_new;
      g_Flu_Array_New[bx][MOMZ][g_idx] = pz_new;

//    for the splitting method, we ensure that the internal energy is unchanged
      g_Flu_Array_New[bx][ENGY][g_idx] = Eint + _rho2*( SQR(px_new) + SQR(py_new) + SQR(pz_new) );
#     endif // #ifdef UNSPLIT_GRAVITY ... else ...


//    update target cell indices
      s_idx_new += NSlice*SQR(GRA_NXT);
#     ifdef UNSPLIT_GRAVITY
      s_idx_old += NSlice*SQR(USG_NXT_G);
#     endif
      g_idx     += NSlice*SQR(PS1);
      z         += NSlice*dh;
   } // for (uint Slice=tz; Slice<PS1; Slice+=NSlice)

} // FUNCTION : CUPOT_HydroGravitySolver



#endif // #if ( defined GPU  &&  MODEL == HYDRO  &&  defined GRAVITY )
