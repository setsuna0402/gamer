#include "GAMER.h"
#include "TestProb.h"


//KH 2023/12/25:
// particle attribute of radiation source 
//========================================================================================
// Bolometric photon number luminosity
static FieldIdx_t Idx_Luminosity   = Idx_Undefined;    // in [#/s] * TimeUnits/LengthUnits^3 (in ENZO)
static FieldIdx_t Idx_CreationTime = Idx_Undefined;    // When Source is formed in code units
static FieldIdx_t Idx_LifeTime     = Idx_Undefined;    // LifeTime of source in code units
static FieldIdx_t Idx_Energy       = Idx_Undefined;    // Photon's Energy[eV] 
// Currently, one RS has one frequency only, to avoid modifying GAMER particle data structure
// static uint N_EnergyBins;       // the number of energy bins
// static FieldIdx_t *Energy_Bins;            // Energy bins [eV] size = N_EnergyBins
// static FieldIdx_t *SED;               // fractional Spectral energy distribution size = N_EnergyBins
static FieldIdx_t Idx_Dir_X        = Idx_Undefined;    // Normalised direction for the north cone of beamed RS 
static FieldIdx_t Idx_Dir_Y        = Idx_Undefined;    // Normalised direction for the north cone of beamed RS 
static FieldIdx_t Idx_Dir_Z        = Idx_Undefined;    // Normalised direction for the north cone of beamed RS
// Currently, GAMER's particle attribute only support FieldIdx_t.
// Change FieldIdx_t Type back to int Type if GAMER supports it. 
// RS->Type: 0:isotropic 1:symmetric beam (north and south pole) 2:non-symmetric beam (north pole)
static FieldIdx_t Idx_Type         = Idx_Undefined;    // Type allows for beaming etc.
static FieldIdx_t Idx_OpenAngle    = Idx_Undefined;    // Opening angle of radiation source in degree
// preserve these valueables for future application
// static FieldIdx_t Idx_RampTime = Idx_Undefined;     // Time for the source to reach full luminosity

//========================================================================================

//KH 2023/12/25:
// particle attribute of radiation source 
//========================================================================================
static FieldIdx_t  Idx_N_Photons_initial    = Idx_Undefined;     // the number of photons created by the RS
static FieldIdx_t  Idx_N_Photons            = Idx_Undefined;     // the actual number of photons in package
// Currently, GAMER's particle attribute only support FieldIdx_t.
// Change FieldIdx_t Type back to int Type if GAMER supports it. 
static FieldIdx_t  Idx_PP_Type              = Idx_Undefined;     // 0 = HI, 1=HeI, 2=HeII, 3=Xray
static FieldIdx_t  Idx_Operation_Type       = Idx_Undefined;     // 0=alive, 1=extinct, 2=sendtoother 3=outofdomain, 4=delete
// static FieldIdx_t  Idx_Energy               = Idx_Undefined;     // Energy of photons in this package [eV]
static FieldIdx_t  Idx_EmissionTimeInterval = Idx_Undefined;     // Duration over which package was emitted
static FieldIdx_t  Idx_BirthTime            = Idx_Undefined;     // Time when package was emitted by RS
static FieldIdx_t  Idx_Radius               = Idx_Undefined;     // Distance travelled
// Currently, GAMER's particle attribute only support FieldIdx_t.
// Change FieldIdx_t HEALPIX_level and HEALPIX_ipix back to int64_t if GAMER supports it. 
static FieldIdx_t  Idx_HEALPIX_level        = Idx_Undefined;     // level in HEALPIX terminology
static FieldIdx_t  Idx_HEALPIX_ipix         = Idx_Undefined;     // pixel in HEALPIX terminology
static FieldIdx_t  Idx_SourceCreationTime   = Idx_Undefined;     // When Radiation source is formed in code units
// Positions are given by OriginalPostion, radius, healpix level and healpix ipix
// static FieldIdx_t  OriginalPosition[3] ;    // Original Position where package was emitted
static FieldIdx_t  Idx_Ori_Pos_X            = Idx_Undefined ;    // Original Position (code unit) where package was emitted
static FieldIdx_t  Idx_Ori_Pos_Y            = Idx_Undefined ;    // Original Position (code unit) where package was emitted
static FieldIdx_t  Idx_Ori_Pos_Z            = Idx_Undefined ;    // Original Position (code unit) where package was emitted
// I don't know whether we need to keep these two parameters in GAMER. They are used to trace PP in my code.
// static FieldIdx_t  Idx_CurrentTime          = Idx_Undefined;     // Current Time
// static FieldIdx_t  Idx_TimeatOneRTTimeStep  = Idx_Undefined;     // how long this PP has traveled at a single RT_timestep
// preseved for future application
// static FieldIdx_t  Idx_SourcePositionDiff = Idx_Undefined;    // Radius at which it was radiated (0 = pt src)

//========================================================================================


// problem-specific global variables
// =======================================================================================
/*
static double ParTest_Dens_Bg;        // background mass density
static double ParTest_Pres_Bg;        // background pressure
static double ParTest_Ang_Freq;       // gas angular frequency
       int    ParTest_NPar[3];        // particles on a side
       double ParTest_Point_Mass;     // the mass of the active particles
       double ParTest_Par_Sep;        // the separation between the active particles
       bool   ParTest_Use_Tracers;    // whether or not to include tracers
       bool   ParTest_Use_Massive;    // whether or not to include massive particles
*/
// =======================================================================================

// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_ParticleTest( const long NPar_ThisRank, const long NPar_AllRank,
                                       real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                       real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                       real *ParType, real *AllAttribute[PAR_NATT_TOTAL] );
#endif

/*
bool Flag_ParticleTest( const int i, const int j, const int k, const int lv,
                        const int PID, const double *Threshold );
*/

//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_RT_Uniform_Rho
// Description :  Add the problem-specific particle attributes
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddParticleField() for each of the problem-specific particle attribute:
//                   --> Attribute label sent to AddParticleField() will be used as the output name of the attribute
//                   --> Attribute index returned by AddParticleField() can be used to access the particle attribute data
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
#if defined(PARTICLE) && defined (RADIATIVE_TRANSER)
void AddNewParticleAttribute_RT_Uniform_Rho()
{
   // RS attributes
   if ( Idx_Luminosity   == Idx_Undefined )
      Idx_Luminosity   = AddParticleAttribute( "Luminosity"   );
   if ( Idx_CreationTime == Idx_Undefined )
      Idx_CreationTime = AddParticleAttribute( "CreationTime" );
   if ( Idx_LifeTime     == Idx_Undefined )
      Idx_LifeTime     = AddParticleAttribute( "LifeTime"     );
   if ( Idx_Energy       == Idx_Undefined )
      Idx_Energy       = AddParticleAttribute( "Energy"       );
   if ( Idx_Dir_X        == Idx_Undefined )
      Idx_Dir_X        = AddParticleAttribute( "Dir_X"        );
   if ( Idx_Dir_Y        == Idx_Undefined )
      Idx_Dir_Y        = AddParticleAttribute( "Dir_Y"        );
   if ( Idx_Dir_Z        == Idx_Undefined )
      Idx_Dir_Z        = AddParticleAttribute( "Dir_Z"        );
   if ( Idx_Type         == Idx_Undefined )
      Idx_Type         = AddParticleAttribute( "Type"         );
   if ( Idx_OpenAngle    == Idx_Undefined )
      Idx_OpenAngle    = AddParticleAttribute( "OpenAngle"    );
   // if ( Idx_RampTime     == Idx_Undefined )
   //    Idx_RampTime     = AddParticleAttribute( "RampTime"     );
   // PP attributes
   if ( Idx_N_Photons_initial    == Idx_Undefined )
      Idx_N_Photons_initial    = AddParticleAttribute( "N_Photons_initial"    );
   if ( Idx_N_Photons            == Idx_Undefined )
      Idx_N_Photons            = AddParticleAttribute( "N_Photons"            );
   if ( Idx_PP_Type              == Idx_Undefined )
      Idx_PP_Type              = AddParticleAttribute( "PP_Type"              );
   if ( Idx_Operation_Type       == Idx_Undefined )
      Idx_Operation_Type       = AddParticleAttribute( "Operation_Type"       );
   // if ( Idx_Energy               == Idx_Undefined )
   //    Idx_Energy               = AddParticleAttribute( "Energy"               );
   if ( Idx_EmissionTimeInterval == Idx_Undefined )
      Idx_EmissionTimeInterval = AddParticleAttribute( "EmissionTimeInterval" );
   if ( Idx_BirthTime            == Idx_Undefined )
      Idx_BirthTime            = AddParticleAttribute( "BirthTime"            );
   if ( Idx_Radius               == Idx_Undefined )
      Idx_Radius               = AddParticleAttribute( "Radius"               );
   if ( Idx_HEALPIX_level        == Idx_Undefined )
      Idx_HEALPIX_level        = AddParticleAttribute( "HEALPIX_level"        );
   if ( Idx_HEALPIX_ipix         == Idx_Undefined )
      Idx_HEALPIX_ipix         = AddParticleAttribute( "HEALPIX_ipix"         );
   if ( Idx_SourceCreationTime   == Idx_Undefined )
      Idx_SourceCreationTime   = AddParticleAttribute( "SourceCreationTime"   );
   if ( Idx_Ori_Pos_X            == Idx_Undefined )
      Idx_Ori_Pos_X            = AddParticleAttribute( "Ori_Pos_X"            );
   if ( Idx_Ori_Pos_Y            == Idx_Undefined )
      Idx_Ori_Pos_Y            = AddParticleAttribute( "Ori_Pos_Y"            );
   if ( Idx_Ori_Pos_Z            == Idx_Undefined )
      Idx_Ori_Pos_Z            = AddParticleAttribute( "Ori_Pos_Z"            );
/*
   if ( Idx_CurrentTime          == Idx_Undefined )
      Idx_CurrentTime          = AddParticleAttribute( "CurrentTime"          );
   if ( Idx_TimeatOneRTTimeStep  == Idx_Undefined )
      Idx_TimeatOneRTTimeStep  = AddParticleAttribute( "TimeatOneRTTimeStep"  );
   if ( Idx_SourcePositionDiff   == Idx_Undefined )
      Idx_SourcePositionDiff   = AddParticleAttribute( "SourcePositionDiff"   );
*/

} // FUNCTION : AddNewParticleAttribute_RT_Uniform_Rho
#endif // #if defined(PARTICLE) && defined (RADIATIVE_TRANSER)

//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

   if ( !OPT__FREEZE_FLUID )
      Aux_Error( ERROR_INFO, "OPT__FREEZE_FLUID must be enabled !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{
/*
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE_ADDRESS,       DEFAULT,      MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "ParTest_Dens_Bg",     &ParTest_Dens_Bg,        1.0e-2,       Eps_double,       NoMax_double      );
   ReadPara->Add( "ParTest_Pres_Bg",     &ParTest_Pres_Bg,        1.0e-2,       Eps_double,       NoMax_double      );
   ReadPara->Add( "ParTest_Ang_Freq",    &ParTest_Ang_Freq,       0.00051668,   Eps_double,       1.0e-3            );
   ReadPara->Add( "ParTest_NParX",       &ParTest_NPar[0],        32,           2,                128               );
   ReadPara->Add( "ParTest_NParY",       &ParTest_NPar[1],        32,           2,                128               );
   ReadPara->Add( "ParTest_NParZ",       &ParTest_NPar[2],        32,           2,                128               );
   ReadPara->Add( "ParTest_Par_Sep",     &ParTest_Par_Sep,        0.5,          Eps_double,       NoMax_double      );
   ReadPara->Add( "ParTest_Point_Mass",  &ParTest_Point_Mass,     1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "ParTest_Use_Tracers", &ParTest_Use_Tracers,    true,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "ParTest_Use_Massive", &ParTest_Use_Massive,    true,         Useless_bool,     Useless_bool      );

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) check the runtime parameters
   if ( !ParTest_Use_Tracers  &&  !ParTest_Use_Massive )
      Aux_Error( ERROR_INFO,
                 "either ParTest_Use_Tracer, ParTest_Use_Massive, or both must be true !!\n" );

#  ifndef TRACER
   if ( ParTest_Use_Tracers )    Aux_Error( ERROR_INFO, "must enable TRACER for ParTest_Use_Tracers !!\n" );
#  endif

#  ifndef GRAVITY
   if ( ParTest_Use_Massive )    Aux_Error( ERROR_INFO, "must enable GRAVITY for ParTest_Use_Massive !!\n" );
#  endif


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 2*M_PI/ParTest_Ang_Freq;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }

// overwrite the total number of particles
#  ifdef PARTICLE
   amr->Par->NPar_Active_AllRank = 0;
   if ( ParTest_Use_Massive )    amr->Par->NPar_Active_AllRank += 2;
   if ( ParTest_Use_Tracers )    amr->Par->NPar_Active_AllRank += ParTest_NPar[0]*ParTest_NPar[1]*ParTest_NPar[2];
   PRINT_WARNING( "PAR_NPAR", amr->Par->NPar_Active_AllRank, FORMAT_LONG );
#  endif


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID            = %d\n",     TESTPROB_ID         );
      Aux_Message( stdout, "  background mass density    = %13.7e\n", ParTest_Dens_Bg     );
      Aux_Message( stdout, "  background pressure        = %13.7e\n", ParTest_Pres_Bg     );
      Aux_Message( stdout, "  angular frequency          = %13.7e\n", ParTest_Ang_Freq    );
      Aux_Message( stdout, "  number of particles (x)    = %d\n",     ParTest_NPar[0]     );
      Aux_Message( stdout, "                      (y)    = %d\n",     ParTest_NPar[1]     );
      Aux_Message( stdout, "                      (z)    = %d\n",     ParTest_NPar[2]     );
      Aux_Message( stdout, "  active particle separation = %13.7e\n", ParTest_Par_Sep     );
      Aux_Message( stdout, "  active particle mass       = %13.7e\n", ParTest_Point_Mass  );
      Aux_Message( stdout, "  include tracer particles   = %d\n",     ParTest_Use_Tracers );
      Aux_Message( stdout, "  include massive particles  = %d\n",     ParTest_Use_Massive );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );
*/
} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

/*
   const double dr[2]     = { x - 0.5*amr->BoxSize[0], y - 0.5*amr->BoxSize[1] };
   const double Radius    = sqrt( dr[0]*dr[0] + dr[1]*dr[1] );

   const double Velocity = ParTest_Ang_Freq*Radius;

   const double Cos_theta = dr[0]/Radius;
   const double Sin_theta = dr[1]/Radius;

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   Dens = ParTest_Dens_Bg;
   Pres = ParTest_Pres_Bg;
   MomX = -ParTest_Dens_Bg*Velocity*Sin_theta;
   MomY = ParTest_Dens_Bg*Velocity*Cos_theta;
   MomZ = 0.0;

   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );    // assuming EoS requires no passive scalars

   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );      // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;
*/
} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_RT_Uniform_Rho
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_RT_Uniform_Rho()
{
/*
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
   Output_User_Ptr               = NULL;
   Flag_User_Ptr                 = Flag_ParticleTest;
   Mis_GetTimeStep_User_Ptr      = NULL;
   Aux_Record_User_Ptr           = NULL;
   BC_User_Ptr                   = NULL;
   Flu_ResetByUser_Func_Ptr      = NULL;
   End_User_Ptr                  = NULL;
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr       = Par_Init_ByFunction_ParticleTest;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
*/
} // FUNCTION : Init_TestProb_RT_Uniform_Rho


/*
bool Flag_ParticleTest( const int i, const int j, const int k, const int lv,
                        const int PID, const double *Threshold )
{

   // Refine a rectanglar solid region in the center just to test behavior in
   // non-uniform regions

   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   const double Center[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double dr[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };

   bool Flag = (FABS(dr[0]) < Threshold[0]) && (FABS(dr[1]) < 2.0*Threshold[0]) && (FABS(dr[2]) < Threshold[0]);

   return Flag;

}
*/

