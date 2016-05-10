#ifndef __CUAPI_H__
#define __CUAPI_H__



#include "Prototype.h"


// CUDA error check 
#define CUDA_CHECK_ERROR( Call )   CUDA_Check_Error( Call, __FILE__, __LINE__, __FUNCTION__ )

inline void CUDA_Check_Error( cudaError Return, const char *File, const int Line, const char *Func )
{
   if ( Return != cudaSuccess ) 
      Aux_Error( ERROR_INFO, "CUDA ERROR : %s !!\n", cudaGetErrorString( Return ) );
}



#endif // #ifndef __CUAPI_H__
