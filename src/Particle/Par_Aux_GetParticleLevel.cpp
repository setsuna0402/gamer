#include "Copyright.h"
#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Aux_GetParticleLevel
// Description :  Record the number of active particles at each level 
//
// Note        :  1. Output filename is "Record__ParLevel"
//                2. It will also check whether the sum of the numbers of active particles
//                   at all levels == Par->NPar_Active
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Par_Aux_GetParticleLevel()
{

   const char FileName[] = "Record__ParLevel";
   static bool FirstTime = true;

// header
   if ( FirstTime )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

      FirstTime = false;

      FILE *File = fopen( FileName, "a" );

      fprintf( File, "#%19s %13s %11s", "Time", "Step", "NActive" );

      for (int lv=0; lv<NLEVEL; lv++)  fprintf( File, " %15s %2d ", "Level", lv );

      fprintf( File, "\n" );

      fclose( File );
   }


// record data
   FILE *File = fopen( FileName, "a" );

   fprintf( File, "%20.13e %13ld %11ld", Time[0], Step, amr->Par->NPar_Active );

   for (int lv=0; lv<NLEVEL; lv++)
      fprintf( File, " %10ld(%6.2lf%%)",
               amr->Par->NPar_Lv[lv], 100.0*(double)amr->Par->NPar_Lv[lv]/(double)amr->Par->NPar_Active );

   fprintf( File, "\n" );

   fclose( File );


// check if sum = amr->Par->NPar_Active
   long sum=0;

   for (int lv=0; lv<NLEVEL; lv++)  sum += amr->Par->NPar_Lv[lv];

   if ( sum != amr->Par->NPar_Active )    Aux_Error( ERROR_INFO, "Sum of active particles (%ld) != expect (%ld) !!\n",
                                                     sum, amr->Par->NPar_Active );

} // FUNCTION : Par_Aux_GetParticleLevel



#endif // #ifdef PARTICLE
