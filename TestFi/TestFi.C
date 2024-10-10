/*--------------------------------------------------------------------------*/
/*-------------------------- File TestFi.C ---------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- TestFi is a simple example of a "concrete" class which implements    --*/
/*-- the interface defined by the abstract base class FiOracle.           --*/
/*--                                                                      --*/
/*--                         Antonio Frangioni                            --*/
/*--                    Dipartimento di Informatica                       --*/
/*--                        Universita' di Pisa                           --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------.-------- INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "TestFi.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------- CONSTANTS ------------------------------------*/
/*--------------------------------------------------------------------------*/

static cIndex InINF = Inf< Index >();

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

TestFi::TestFi( cIndex NV , cLMNum L0
		//!! these two parameters are only needed for the sake of
		//!! this example: remove them and replace with your own
		)
        :
	FiOracle()
        // constructors of other ancestor classes, if any
{
 // initializations, if needed, possibly making use of the extra parameters;
 // allocation of TestFi-specific data structures is usually performed here

 //!! the following code is needed only for the sake of this example
 //!! remove it and replace with your own

 NumVar = NV;
 Cntr = L0;

 //!! end of example-specific code
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

HpNum TestFi::Fi( cIndex wFi )
{
 if( wFi ) {
  //!! the following code is needed only for the sake of this example:
  //!! compute
  //!!
  //!!  Fi( Lambda ) = 1/2 * sum{ i = 1 .. NV } ( Lambda[ i ] - Cntr )^2 + 1
  //!!
  //!! remove it and replace with your own

  HpNum temp = 0;

  if( LamBase )
   for( Index i = 0 , j = 0 ; i < NumVar ; i++ ) {
    HpNum SGi = - Cntr;
    if( i == LamBase[ j ] )
     SGi += Lambda[ j++ ];

    temp += SGi * SGi;
    }
  else
   for( Index i = 0 ; i < NumVar ; i++ ) {
    cHpNum SGi = ( Lambda[ i ] - Cntr );
    temp += SGi * SGi;
    }

  return( temp / 2 + 1 );

  //!! end of example-specific code
  }
 else
  return( 0 );

 }  // end( TestFi::Fi )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

Index TestFi::GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name ,
		     cIndex strt , Index stp )
{
 if( Name < MaxName )
  throw( NDOException( "GetGi: past information is not recorded" ) );

 if( strt || ( stp <= NumVar ) )
  throw( NDOException( "GetGi: slicing the subgradient is not supported" ) );

 if( Name == MaxName ) {  // the 0-th component
  SGBse = &InINF;         // ... is an all-0 vector
  return( 0 );
  }

 //!! the following code is needed only for the sake of this example:
 //!! compute
 //!!
 //!!   Gi( Lambda )[ i ] = Lambda[ i ] - Cntr
 //!!
 //!! remove it and replace with your own

 if( LamBase )
  for( Index i = 0 , j = 0 ; i < NumVar ; i++ )
   if( i == LamBase[ j ] )
    SubG[ i ] = ( Lambda[ j++ ] - Cntr );
   else
    SubG[ i ] = - Cntr;
 else
  for( Index i = 0 ; i < NumVar ; i++ )
   SubG[ i ] = ( Lambda[ i ] - Cntr );

 //!! end of example-specific code

 LHasChgd = false;  // important: only one subgradient for each Lambda[]
 SGBse = 0   ;      // a "dense" subgradient ...
 return( NumVar );  // ... has NumVar nonzero entries (and if some are
                    // actually zero, we don't care)
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

TestFi::~TestFi()
{
 // deallocation of TestFi-specific data structures (belonging to TestFi but
 // *not* to the base classes, like FiOracle) is typically performed here
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------ End File TestFi.C -------------------------------*/
/*--------------------------------------------------------------------------*/
