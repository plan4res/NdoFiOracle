/*--------------------------------------------------------------------------*/
/*--------------------------- File test.cpp --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Simple main() for testing the solvers of the NODSolver / FiOracle using
 * the TestFi oracle.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#define HAVE_CutPlane 0
// if nonzero, the CutPlane NDOSolver is loaded

#define HAVE_OsiMPSolver 1
// if nonzero, the OsiMPSolver MPSolver is loaded

#if HAVE_CutPlane || HAVE_OsiMPSolver
 #define HAVE_OsiSolverInterface 1
#else
 #define HAVE_OsiSolverInterface 0
#endif
// OsiSolverInterface only need to be loaded if either HAVE_CutPlane or
// HAVE_OsiMPSolver is nonzero

#if HAVE_OsiSolverInterface
 #define HAVE_OsiCpxSolverInterface 1
 // if nonzero, the OsiCpxSolverInterface is loaded; this makes only sense if
 // OsiSolverInterface is used
#endif

#define PRINT_SOLUTION 1
// if nonzero, the "optimal" point is printed in its full glory

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

// the FiOracle- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "TestFi.h"

// the CutPlane solver - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#if HAVE_CutPlane
 #include "CutPlane.h"
#endif

// the Bundle solver - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#include "Bundle.h"
#include "QPPnltMP.h"

#if HAVE_OsiMPSolver
  #include "OSIMPSolver.h"
#endif

#if HAVE_OsiSolverInterface
 #include "OsiClpSolverInterface.hpp"
 // OsiClpSolverInterface is "free" and therefore is included by default

 #if HAVE_OsiCpxSolverInterface
  #include "OsiCpxSolverInterface.hpp"
  #include "ilcplex/cplex.h"
 #endif
#endif

// the SubGrad solver- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include "SubGrad.h"

// Deflection headers

#include "Volume.h"
#include "PrimalDual.h"

// Stepsize headers

#include "ColorTV.h"
#include "FumeroTV.h"
#include "Polyak.h"

// other headers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#include <fstream>
#include <sstream>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

template< class T >
static inline void str2val( const char* const str , T &sthg )
{
 istringstream( str ) >> sthg;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/
/* The configuration of the NDOSolver is completely described by the values
 * found in the parameter file. The format is:
 *
 * - which NDOSolver is used [integer, default 1]:
 *   0 = CutPlane [if available, otherwise error is reported]
 *   1 = Bundle
 *   2 = SubGrad
 *
 * - verbosity of the NDOSolver [integer, defaut 0]
 *
 * After that, the parameter file depends on the previous choices.
 *
 * - - - - - For CutPlane:
 *
 * - which OsiXXXSolverInterface is used [integer, defaut 0]:
 *   0 = OsiClpSolverInterface
 *   1 = OsiCpxSolverInterface [if available, otherwise error is reported]
 *
 * after which the configuration file is passed to the CutPlane constructor
 * (which in turn passes it first to the NDOSolver constructor), see the
 * comments there for details.
 *
 * - - - - - For Bundle:
 *
 * - which MPSolver is used [integer, defaut 0]:
 *   0 = QPPenaltyMP
 *   1 = OsiMPSolver with boxstep stabilization [if available, else ignored]
 *   2 = OsiMPSolver with proximal stabilization [if available, else ignored]
 *
 * - which OsiXXXSolverInterface is used [integer, defaut 0] if OsiMPSolver
 *   is used, ignored otherwise:
 *   0 = OsiClpSolverInterface
 *   1 = OsiCpxSolverInterface [if available, otherwise ignored]
 *
 * after which the configuration file is passed, in this order:
 * = to the Bundle constructor (which in turn passes it first to the NDOSolver
 *   constructor)
 * = to the specific MPSolver constructor
 * see the comments there for details.
 *
 * - - - - - For SubGrad:
 *
 * - which Deflection is used [integer, defaut 0]:
 *   0 = none
 *   1 = Volume
 *   2 = PrimalDual
 *
 * - which Stepsize is used [integer, defaut 0]:
 *   0 = ColorTV
 *   1 = FumeroTV
 *   2 = Polyak
 *   3 = PrimalDual (note that this should be used if and only if PrimalDual
 *                   is used as Deflection)
 *
 * after which the configuration file is passed, in this order:
 * = to the SubGrad constructor (which in turn passes it first to the
 *   NDOSolver constructor)
 * = to the specific Deflection constructor (if used)
 * = to the specific Stepsize constructor (if used)
 * see the comments there for details.
 */

int main( int argc , char **argv )
{
 // read the command-line parameters- - - - - - - - - - - - - - - - - - - - -

 if( argc < 4 ) {
  cerr << "Usage: " << argv[ 0 ]
       << " < # variables > < center > < parfile name >" << endl;

  return( 1 );
  }

 Index NV;
 str2val( argv[ 1 ] , NV );

 LMNum Cntr;
 str2val( argv[ 2 ] , Cntr );

 ifstream ParFile( argv[ 3 ] );
 if( ! ParFile.is_open() ) {
  cerr << "Cannot open parameters file """ << argv[ 3 ] << """" << endl;
  return( 1 );
  }

 int WNdo;
 DfltdSfInpt( &ParFile , WNdo , int( 0 ) );
 if( ( WNdo < 0 ) || ( WNdo > 2 ) ) {
  cerr << "NDO parameter " << WNdo << " invalid" << endl;
  return( 1 );
  }

 int lvl;
 DfltdSfInpt( &ParFile , lvl , int( 0 ) );

 int optn1;
 DfltdSfInpt( &ParFile , optn1 , int( 0 ) );

 // enter the try-block - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 try {
  // construct the FiOracle object- - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FiOracle *Fi = new TestFi( NV , Cntr );

  // construct the NDOSolver object - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  NDOSolver *NDOslvr;
  MPSolver *MPslvr = nullptr;
  #if HAVE_OsiSolverInterface
   OsiSolverInterface *OSIslvr = nullptr;
  #endif

  switch( WNdo ) {
   case( 0 ): {  // CutPlane - - - - - - - - - - - - - - - - - - - - - - - - -

    #if HAVE_CutPlane
     CutPlane *ctpln = new CutPlane( &ParFile );

     #if HAVE_OsiCpxSolverInterface
      if( optn1 ) {
       // need to use the methods of OsiCpxSolverInterface to retrieve the
       // Cplex environment to set some parameters because this is not easy
       // to do via OsiSolverInterface
       OsiCpxSolverInterface *tos = new OsiCpxSolverInterface();

       CPXENVptr env = tos->getEnvironmentPtr();

       // try to ensure that the solver stays silent
       CPXsetintparam( env , CPX_PARAM_SCRIND , CPX_OFF );
       CPXsetintparam( env , CPX_PARAM_MIPDISPLAY , CPX_OFF );
       CPXsetintparam( env , CPX_PARAM_PARAMDISPLAY , CPX_OFF );
       CPXsetintparam( env , CPXPARAM_ScreenOutput , CPX_OFF );
       CPXsetintparam( env , CPXPARAM_Barrier_Display , 0 );
       CPXsetintparam( env , CPXPARAM_Simplex_Display , 0 );
       CPXsetintparam( env , CPXPARAM_Sifting_Display , 0 );
       CPXsetintparam( env , CPXPARAM_Network_Display , 0 );
       CPXsetintparam( env , CPXPARAM_ParamDisplay , CPX_OFF );
 
       // force it to run with only one thread
       CPXsetintparam( env , CPX_PARAM_THREADS , 1 );
       OSIslvr = tos;
       }
      else
     #endif
      {
       OsiClpSolverInterface *tos = new OsiClpSolverInterface();
       tos->getModelPtr()->setLogLevel( 0 );
       OSIslvr = tos;
       }

     OSIslvr->setHintParam( OsiDoReducePrint );

     ctpln->SetOsiSolver( OSIslvr );
     //    ctpln->SetRadius( 10 * Cntr );
     NDOslvr = ctpln;
     break;
    #else
     cout << "CutPlane solver not available" << endl;
     return( 1 );
    #endif
    }
   case( 1 ): {  // Bundle - - - - - - - - - - - - - - - - - - - - - - - - - -

    int optn2;
    DfltdSfInpt( &ParFile , optn2 , int( 0 ) );

    Bundle *bndl = new Bundle( &ParFile );

    if( optn1 ) {
     #if HAVE_OsiMPSolver
      OSIMPSolver *osimp = new OSIMPSolver();

      if( optn2 ) {
       #if HAVE_OsiCpxSolverInterface
        // need to use the methods of OsiCpxSolverInterface to retrieve the
        // Cplex environment to set some parameters because this is not easy
        // to do via OsiSolverInterface
        OsiCpxSolverInterface *tos = new OsiCpxSolverInterface();

        CPXENVptr env = tos->getEnvironmentPtr();

	// try to ensure that the solver stays silent
	CPXsetintparam( env , CPX_PARAM_SCRIND , CPX_OFF );
	CPXsetintparam( env , CPX_PARAM_MIPDISPLAY , CPX_OFF );
	CPXsetintparam( env , CPX_PARAM_PARAMDISPLAY , CPX_OFF );
	CPXsetintparam( env , CPXPARAM_ScreenOutput , CPX_OFF );
	CPXsetintparam( env , CPXPARAM_Barrier_Display , 0 );
	CPXsetintparam( env , CPXPARAM_Simplex_Display , 0 );
	CPXsetintparam( env , CPXPARAM_Sifting_Display , 0 );
	CPXsetintparam( env , CPXPARAM_Network_Display , 0 );
	CPXsetintparam( env , CPXPARAM_ParamDisplay , CPX_OFF );

        // force it to run with only one thread
        CPXsetintparam( env , CPX_PARAM_THREADS , 1 );
        OSIslvr = tos;
       #else
        cout << "OsiCpxSolverInterface not available" << endl;
        return( 1 );
       #endif
       }
      else {
       OsiClpSolverInterface *tos = new OsiClpSolverInterface();
       tos->getModelPtr()->setLogLevel( 0 );
       OSIslvr = tos;
       }

      OSIslvr->setHintParam( OsiDoReducePrint );

      osimp->SetOsi( OSIslvr );
      osimp->SetStabType( optn1 == 1 ? OSIMPSolver::boxstep
			             : OSIMPSolver::quadratic );
      MPslvr = osimp;
     #else
      cout << "OSIMPSolver not available" << endl;
      return( 1 );
     #endif
     }
    else
     MPslvr = new QPPenaltyMP( &ParFile );

    bndl->SetMPSolver( MPslvr);
    NDOslvr = bndl;
    break;
    }
   case( 2 ): {  // SubGrad- - - - - - - - - - - - - - - - - - - - - - - - - -

    int optn2;
    DfltdSfInpt( &ParFile , optn2 , int( 0 ) );

    SubGrad *sbgrd = new SubGrad( &ParFile );

    Deflection *dflctn = nullptr;
    switch( optn1 ) {
     case( 0 ): break;
     case( 1 ): dflctn = new Volume( sbgrd , &ParFile ); break;
     case( 2 ): dflctn = new PrimalDual( sbgrd , &ParFile ); break;
     default:
      cerr << "Invalid Deflection " << optn1 << endl; return( 1 );
     }

    Stepsize *stpsz;
    switch( optn2 ) {
     case( 0 ): stpsz = new ColorTV( sbgrd , &ParFile ); break;
     case( 1 ): stpsz = new FumeroTV( sbgrd , &ParFile ); break;
     case( 2 ): stpsz = new Polyak( sbgrd , &ParFile ); break;
     case( 3 ): stpsz = dynamic_cast< PrimalDual * >( dflctn );
      if( stpsz )
       break;
      else {
       cerr << "PrimalDual Stepsize requires PrimalDual Deflection" << endl;
       return( 1 );
       }
     default:
      cerr << "Invalid Stepsize " << optn2 << endl; return( 1 );
     }

    sbgrd->SetStepsize( stpsz );
    sbgrd->SetDeflection( dflctn );
    NDOslvr = sbgrd;
    }
   }  // end( switch( WNdo ) ): construct the NDOSolver - - - - - - - - - - -
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // close the parameters file  - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ParFile.close();

  // pass the FiOracle to the NDOSolver - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  NDOslvr->SetFiOracle( Fi );

  // set the verbosity of log - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( lvl )
   NDOslvr->SetNDOLog( &clog , char( lvl ) );

  // minimize the function- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  NDOslvr->Solve();

  // print results- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << endl << "Fi* = " << NDOslvr->ReadBestFiVal();
  if( NDOslvr->IsOptimal() )
   cout << " (optimal value)" << endl;
  else
   cout << " (not provably optimal)" << endl;

  #if PRINT_SOLUTION
   cout << "Lambda* =" << endl;
   Index D;
   cIndex_Set I;
   cLMRow L = NDOslvr->ReadSol( I , D );

   if( I )
    for( Index i = 0 ; ( i = *(I++) ) < Inf< Index >() ; )
     cout << "[" << i << "]\t" << *(L++) << endl;
   else
    for( Index i = 0 ; i < NDOslvr->GetNumVar() ; i++ )
     cout << "[" << i << "]\t" << *(L++) << endl;
  #endif

  // delete the objects - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( lvl )
   NDOslvr->SetNDOLog( nullptr , 0 );

  delete( NDOslvr );

  delete( MPslvr );

  #if HAVE_OsiSolverInterface
   delete( OSIslvr );
  #endif

  delete( Fi );

  }  // end( try-block )

 // managing exceptions - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 catch( exception &e ) {
  cerr << e.what() << endl;
  return( 1 );
  }
 catch(...) {
  cerr << "Error: unknown exception thrown" << endl;
  return( 1 );
  }

 // the end - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( 0 );

 }  // end( main )

/*--------------------------------------------------------------------------*/
/*------------------------ End File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
