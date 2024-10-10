/*--------------------------------------------------------------------------*/
/*----------------------------- File CutPlane.C ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the CutPlane class, which implements the NDOSolver 
 * interface for NonDifferentiable Optimization Solvers, as described in
 * NDOSlver.h, using the "Cutting Plane" algorithm.                     
 * Uses a generic OsiSolverInterface to "communicate" with the linear
 * programming solver used to solve the (restricted) master problem.
 *                                                                      
 * \author Vitor Barbosa \n
 *         Algoritmi Research Unit \n
 *         Universidade do Minho \n
 *
 * \author Filipe Alvelos \n
 *         Grupo de Optimizacao e Investigacao Operacional \n
 *         Departamento de producao e Sistemas \n
 *         Universidade do Minho \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Filipe Alvelos, Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CutPlane.h"

#include "CoinPackedVectorBase.hpp"

#include "CoinPackedVector.hpp"

#include "OPTvect.h"

#include <algorithm>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_CUT 1
/* If LOG_CUT > 0, the CutPlane class may produce (depending on an input
   parameter) a log of (more or less) important results. */

#if( LOG_CUT )
 #define CUTLOG( y , x ) if( NDOLLvl >= y ) *NDOLog << x
 #define CUTLOG2( y , c , x ) if( ( NDOLLvl >= y ) && c ) *NDOLog << x
#else
 #define CUTLOG( y , x ) 
 #define CUTLOG2( y , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define RMPAlg 2

/* Algorithm to solve the Restricted Master Problems
    1 Use OsiSolverInterface function initialSolve() 
    2 Use OsiSolverInterface function Resolve()
    */

/*--------------------------------------------------------------------------*/
/*---------------------- IMPLEMENTATION OF CutPlane ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

CutPlane::CutPlane( istream *iStrm ) : NDOSolver( iStrm )
{
 // initialize algorithmic parameters - - - - - - - - - - - - - - - - - - - -

 DfltdSfInpt( iStrm , NumNames , Index( 100 ) );
 DfltdSfInpt( iStrm , PurgeRows , HpNum( -1 ) );
 DfltdSfInpt( iStrm , PurgeInvl , Index( 5 ) );
 DfltdSfInpt( iStrm , HeurSubGU , Index( 1 ) );
 DfltdSfInpt( iStrm , HeurSubGF , Index( 1 ) );
 DfltdSfInpt( iStrm , HeurInvl , Index( 30 ) );
 DfltdSfInpt( iStrm , KpPrimals , Index( 0 ) );
 DfltdSfInpt( iStrm , PrintInvl , Index( 50 ) );

 // some initializations- - - - - - - - - - - - - - - - - - - - - - - - - - -

 KpBstL = false;
 BestLambda = 0;
 AllocNumber = 0;
 OsiS = 0;

 }  // end( CutPlane )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void CutPlane::SetOsiSolver( OsiSolverInterface* osi )
{
 if( Oracle )
  throw( NDOException( "CutPlane::SetOsiSolver: FiOracle already present" ) );

 OsiS = osi;
 }  

/*--------------------------------------------------------------------------*/

void CutPlane::SetRadius( LMNum rad )
{
 const double OsiINF = OsiS->getInfinity();

 if( rad == Inf< LMNum >() ) {  // resetting the box- - - - - - - - - - - - -
  for( Index i = 0 ; i < NumVar ; i++ )
   OsiS->setColBounds( NrFi + i ,
		       Oracle->GetUC( i ) ? - OsiINF : double( 0 ) ,
		       Oracle->GetUB( i ) == Inf< LMNum >() ? OsiINF :
		                              double( Oracle->GetUB( i ) ) );
  }
 else {                // setting the box - - - - - - - - - - - - - - - - - -
  if( BaseSize < NumVar ) {
   cLMRow tL = Lambda;
   cIndex_Set tB = Base;
   for( Index i = 0 ; i < NumVar ; i++ ) {
    LMNum Li = 0;
    if( i == *tB ) {
     Li = *(tL++);
     tB++;
     }
    double ubi = Oracle->GetUB( i ) == Inf< LMNum >() ? OsiINF :
                                                      Oracle->GetUB( i );
    ubi = std::min( ubi , double( Li + rad ) );
    double lbi = Oracle->GetUC( i ) ? - OsiINF : 0;
    lbi = std::max( lbi , double( Li - rad ) );
    OsiS->setColBounds( NrFi + i , lbi , ubi );
    }
   }
  else
   for( Index i = 0 ; i < NumVar ; i++ ) {
    double ubi = Oracle->GetUB( i ) == Inf< LMNum >() ? OsiINF :
                                                      Oracle->GetUB( i );
    ubi = std::min( ubi , double( Lambda[ i ] + rad ) );
    double lbi = Oracle->GetUC( i ) ? - OsiINF : 0;
    lbi = std::max( lbi , double( Lambda[ i ] - rad ) );
    OsiS->setColBounds( NrFi + i , lbi , ubi );
    }
  }
 }  // end( CutPlane::SetRadious )

/*--------------------------------------------------------------------------*/

void CutPlane::SetFiOracle( FiOracle *Fi )
{
 if( Oracle )  // changing from a previous oracle - - - - - - - - - - - - - -
 {             // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  MemDealloc();   // deallocate memory
  }

 if( Fi )      // setting a new oracle- - - - - - - - - - - - - - - - - - - -
 {             // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ! OsiS )
  throw( NDOException( "CutPlane::SetFiOracle: OsiSolver not present" ) );

  // "throw" the method of the base class - - - - - - - - - - - - - - - - - -

  NDOSolver::SetFiOracle( Fi );

  // tell the oracle about the NDOSolver and its settings - - - - - - - - - -

  Oracle->SetNDOSolver( this );
  Oracle->SetPrecision( EFnal );

  // initialize GAP

  Gap = Inf< HpNum >();

  // read information about Fi- - - - - - - - - - - - - - - - - - - - - - - -

  MaxNumVar = Oracle->GetMaxNumVar();
  LowerBound = Oracle->GetLowerBound();
  if( LowerBound > -Inf< HpNum >() )
   TrueLB = true;
  else {
   TrueLB = false;
   LowerBound = Oracle->GetMinusInfinity();
   }

  // allocate memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Lambda = new LMNum[ MaxNumVar ];  // current point
  Base = new Index[ MaxNumVar ];    // indices of nonzeroes in Lambda
  Gi = new SgNum[ MaxNumVar ];      // subgradient
  FiLambda = new HpNum[ NrFi + 1 ];
  BestFiLambda = new HpNum[ NrFi + 1 ];
  ValueRMP = new HpNum[ NrFi + 1 ];

  // allocate memory for "names" management

  if( KpPrimals ) {
   FullNames = new Index[ NumNames ];
   SubPNames = new Index[ NumNames ];
   EmptyNames = new Index[ NumNames ];
   FirstEmpty = 0;
   LastFull = -1;
   AllocNumber = 1;
   for( Index i = 0 ; i < NumNames ; i++ )
    EmptyNames[ i ] = i;

   Oracle->SetMaxName( NumNames );   }
  else {
   FullNames = SubPNames = EmptyNames = 0;
   Oracle->SetMaxName( 0 );
   }

  // some initializations - - - - - - - - - - - - - - - - - - - - - - - - - -

  SetLambda( 0 );
  for( Index i = 0 ; i <= NrFi ; i++ ) {
   ValueRMP[ i ] = - Inf< HpNum >();
   BestFiLambda[ i ] = Inf< HpNum >();
   }

  AValueRMP = - Inf< HpNum >();
  ABestFiLambda = Inf< HpNum >();

  Result = kError;
  ParIter = 0;

  // now construct the initial RMP- - - - - - - - - - - - - - - - - - - - - -
  // the RMP is:
  //
  //  min sum_{i=1}^NrFi v[ i ] + g0 * y
  //      sum_{i=1}^NrFi v[ i ] + g0 * y >= LB
  //      v[ i ] >= y * g[ j ] + alpha[ j ]   j \in B( i ) , i = 1, ..., NrFi
  //      l[ h ] <= y[ h ] <= u[ h ]          h = 0, ..., NumVar 
  //
  // where B( i ) is the set of subgradients corresponding to the i-th
  // component of the function, and g0 is the (sub)gradient of the 0-th
  // component. note the first constraint, which ensures boundedness below
  // of the RMP if at least some rough lower bound can be computed

  // objective function sense

  OsiS->setObjSense( 1 );

  // empty CoinPackedVector

  // get the Osi INF value once and for all

  const double OsiINF = OsiS->getInfinity();

  // add "empty" columns for the v[ i ] variables to the OsiS

  for( Index i = 0 ; i < NrFi ; i++ )
   OsiS->addCol( 0 , 0 , 0 , - OsiINF , OsiINF , 1 );

  // add "empty" lambda columns to the OsiS; for now assume g0 = 0, if
  // not this will be changed below

  for( Index i = 0 ; i < NumVar ; i++ )
   OsiS->addCol( 0 , 0 , 0 ,
		 Oracle->GetUC( i ) ? - OsiINF : 0 ,
	         Oracle->GetUB( i ) == Inf< LMNum >() ? OsiINF :
		                                      Oracle->GetUB( i ) ,
		 0 );

  // add global lowerbound constraint

  CoinPackedVector vect;
  for( Index i = 0 ; i < NrFi ; i++ )
   vect.insert( i , 1 );

  // get objective function coeficients from the oracle

  cIndex_Set SGBse;
  cIndex SGBDm = Oracle->GetGi( Gi , SGBse , Oracle->GetMaxName() );

  if( SGBse )
   for( Index i = 0 ; i < SGBDm ; i++ ) {
    vect.insert( SGBse[ i ] + NrFi , Gi[ i ] );
    OsiS->setObjCoeff( SGBse[ i ] + NrFi , Gi[ i ] );
    }
  else
   for( Index i = 0 ; i < NumVar ; i++ ) {
    vect.insert( i + NrFi , Gi[ i ] );
    OsiS->setObjCoeff( i + NrFi , Gi[ i ] );
    }

  OsiS->addRow( vect , LowerBound > -Inf< HpNum >() ? LowerBound : - OsiINF ,
		                                                   OsiINF );

  }  // end( if( Fi ) ) - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MaxNRows = 0;

 }  // end( CutPlane::SetFiOracle )

/*--------------------------------------------------------------------------*/

void CutPlane::SetLambda( cLMRow tLambda )
{
 if( ! Oracle )
  throw( NDOException( "CutPlane::SetLambda: Oracle == 0" ) );

 if( tLambda ) {
  for( Index i = 0 ; i < NumVar ; i++ ) {
   if( ( ! Oracle->GetUC( i ) ) && ( tLambda[ i ] < 0 ) )
    throw( NDOException( "CutPlane::SetLambda: lower bound infeasible" ) );

   if( tLambda[ i ] > Oracle->GetUB( i ) )
    throw( NDOException( "CutPlane::SetLambda: upper bound infeasible" ) );

   Lambda[ i ] = tLambda[ i ];
   }

  Index_Set tB = Sparsify( Lambda , Base , NumVar );
  if( ( BaseSize = tB - Base ) < NumVar )
   *tB = Inf< Index >();
  }
 else {
  for( Index i = 0 ; i < NumVar ; i++ ) {
   if( Oracle->GetUB( i ) < 0 )
    throw( NDOException( "CutPlane::SetLambda: upper bound infeasible" ) );

   Lambda[ i ] = 0;
   }

  *Base = Inf< Index >();
  BaseSize = 0;
  }

 AFiLambda = Inf< HpNum >();
 VectAssign( FiLambda , HpNum( Inf< HpNum >() ) , NrFi + 1 );

 }  // end( SetLambda )

/*--------------------------------------------------------------------------*/

void CutPlane::KeepBestLambda( const bool KBL )
{
 if( KBL && ( ! KpBstL ) )
  BestLambda = new LMNum[ NumVar ]; // best point

 if( ( ! KBL ) && KpBstL ) {
  delete[] BestLambda;
  BestLambda = 0;
  }

 KpBstL = KBL;
 }  
  
/*--------------------------------------------------------------------------*/

void CutPlane::SetNDOLog( ostream *outs, const char lvl )
{
 if( ( NDOLog = outs ) )
  NDOLLvl = lvl;
 else
  NDOLLvl = 0;
 }
  
/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

NDOSolver::NDOStatus CutPlane::Solve( void )
{
 // initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NDOt )
  NDOt->Start();

 Result = kOK;
 SCalls++;        // increment Solve() calls counter

 // initialization step - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! ParIter ) {
  FiAndGi();    // get first subgradient from the Oracle

  if( Result )  // something was wrong in the Fi() computation
   return( Result );

  // write first Fi info to log file

  CUTLOG( 2 , std::endl << "0 :\tFi(): " << AFiLambda << std::endl );

  // check the Lower Bound- - - - - - - - - - - - - - - - - - - - - - - - - -

  if( AFiLambda < LowerBound ) {
   Result = kUnbndd;
   CUTLOG( 0 , std::endl << "\tSTOP: Fi unbounded " << std::endl );
   return( Result );
   }
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for(;;) {
  RowsAdded = 0;  // set rows added in current iteration to 0
  
  // solve the Restricted Master Problem- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( RMPAlg == 1 )
   OsiS->initialSolve();
  #else
   OsiS->resolve();
  #endif

  MaxNRows = std::max( MaxNRows , Index( OsiS->getNumRows() ) );

  if( ! OsiS->isProvenOptimal() ) {
   CUTLOG( 0 , std::endl << "RMP has no optimal solution" << std::endl );
   Result = kError;
   break;
   }

  const double *SolRMP = OsiS->getColSolution();  // get solution values
  AValueRMP = OsiS->getObjValue();                // get objective value
  Gap = ABestFiLambda - AValueRMP;

  HpNum Sum_U = 0;
  for( Index i = 1 ; i <= NrFi ; i++ ) {
   ValueRMP[ i ] = SolRMP[ i - 1 ];
   Sum_U += ValueRMP[ i ];
   }

  ValueRMP[ 0 ] = AValueRMP - Sum_U;

  // set lambda equal to the optimal solution of the RMP

  for( Index i = 0 ; i < NumVar ; i++ )
   Lambda[ i ] = SolRMP[ i + NrFi ];

  Index_Set tB = Sparsify( Lambda , Base , NumVar );
  if( ( BaseSize = tB - Base ) < NumVar )
   *tB = Inf< Index >();

  // check the status of the FiOracle and take the necessary action - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FiOracle::FiStatus fs = Oracle->GetFiStatus();
  if( fs == FiOracle::kFiStop ) {
   CUTLOG( 0 , " ~ FiOracle:STOP" << std::endl );
   Result = kStopped;
   break;
   }

  if( fs == FiOracle::kFiError ) {
   CUTLOG( 0 , " ~ Error in the FiOracle" << std::endl );
   Result = kError;
   break;
   }

  if( fs == FiOracle::kFiChgd ) {
   CUTLOG( 1 , " ~ Fi changed: loop" << std::endl );
   continue;
   }

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( fs != FiOracle::kFiCont )
   if( CutPlane::IsOptimal() )
    break;

  // check for running time - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( MaxTime && NDOt )
   if( NDOt->Read() > MaxTime ) {
    Result = kStopped;
    break;
    }

  ParIter++;  // increment iterations counter

  // delete rows- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index DeletedRows = 0;
  if( ( PurgeRows >= 0 ) && ( ParIter % PurgeInvl == 0 ) )
   DeletedRows = DelRMPRows( std::max( PurgeRows * Gap , 1e-6 ) );

  // calculate Fi( Lambda ) - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FiAndGi();

  Gap = ABestFiLambda - AValueRMP;

  // print some numbers - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( LOG_CUT )
   if( ( NDOLLvl > 2 ) ||
       ( ( NDOLLvl == 2 ) && ( ! ( ParIter % PrintInvl ) ) ) )
    *NDOLog << std::endl << ParIter << ":\tValueRMP: " << AValueRMP
	    << " ~ Fi(): " << AFiLambda << " ~ Gap: " << Gap << std::endl
	    << "\tAddRows: " << RowsAdded << " ~ DelRows: " << DeletedRows
	    << " ~ NumRows: " << OsiS->getNumRows() << std::endl;
  #endif

  if( Result )  // something was wrong in the Fi() computation
   break;

  // re-check for optimality- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( fs != FiOracle::kFiCont )
   if( CutPlane::IsOptimal() )
    break;

  // check that something will change - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( RowsAdded == 0 ) {
   Result = kError;
   CUTLOG( 0 , std::endl << "\tSTOP: no new subgradients added"
	                 << std::endl );
   break;
   }

  // check the Lower Bound- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( AFiLambda < LowerBound ) {
   Result = kUnbndd;
   CUTLOG( 0 , std::endl << "\tSTOP: Fi unbounded " << std::endl );
   break;
   }

 // check iterations count - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( MaxIter && ( ParIter >= MaxIter ) ) {
   Result = kStpIter;
   break;
   }
  }  // end( for( ever ) )

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NDOt )
  NDOt->Stop();

 OutRsts();

 return( Result );

 }  // end( CutPlane::Solve )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

cLMRow CutPlane::ReadSol( cIndex_Set &I , Index &D )
{
 I = BaseSize < NumVar ? Base : 0;
 D = BaseSize;
 return( Lambda );
 }

/*--------------------------------------------------------------------------*/

cLMRow CutPlane::ReadBestSol( cIndex_Set &I , Index &D )
{
 I = 0;
 D = NumVar;
 return( KpBstL ? BestLambda : 0 );
 }

/*--------------------------------------------------------------------------*/

HpNum CutPlane::ReadFiVal( cIndex wFi )
{
 if( wFi == Inf< Index >() )
  return( AFiLambda );
 else
  if( wFi > NrFi )
   return( SumV( FiLambda + 1 , NrFi ) );
  else
   return( FiLambda[ wFi ] );
 }

/*--------------------------------------------------------------------------*/

HpNum CutPlane::ReadBestFiVal( cIndex wFi )
{
 if( wFi == Inf< Index >() )
  return( ABestFiLambda );
 else
  if( wFi > NrFi )
   return( SumV( BestFiLambda + 1 , NrFi ) );
  else
   return( BestFiLambda[ wFi ] );
 }

/*--------------------------------------------------------------------------*/

bool CutPlane::IsOptimal( HpNum eps ) const
{
 if( eps == 0 )
  eps = EpsLin;

 return( ( ABestFiLambda < Inf< HpNum >() ) &&
	 ( Gap <= eps * std::max( ABS( ABestFiLambda ) , HpNum( 1 ) ) ) );
 }

/*--------------------------------------------------------------------------*/

cHpRow CutPlane::ReadMult( cIndex_Set &I , Index &D , cIndex wFi )
{
 if( ! KpPrimals )
  throw( NDOException( "ReadMult Error: past info not recorded." ) );

 if( ( wFi != Inf< Index >() ) && ( wFi > NrFi) )
  throw( NDOException( "ReadMult Error: subproblem dont exists." ) );

  D = 0;
  OsiS->resolve();
  cHpRow duals = OsiS->getRowPrice();
  HpRow res;
  Index_Set vals;
    
  if( wFi == Inf< Index >() ) {
   D = OsiS->getNumRows() - 1;
   res = new HpNum[ D ];
   vals = new Index[ D ];
   for( Index i = 0 ; i < D ; i++ ) {
    res[ i ] = duals[ i + 1 ];
    vals[ i ] = FullNames[ i ];
    }
   }
  else {
   for( Index i = 0 ; i < Index( OsiS->getNumRows() - 1 ) ; i++ )
    if( SubPNames[ i ] == wFi )
     D++;

   res = new HpNum[ D ];
   vals = new Index[ D ];
   Index pos = 0;
   for( Index i = 0 ; i < Index( OsiS->getNumRows() - 1 ) ; i++ ) {
    if( SubPNames[ i ] == wFi ) {
     res[ pos ] = duals[ i  + 1 ];
     vals[ pos ] = FullNames[ i ];
     pos++;
     }
    }
   }

 I = vals;
 return( res );

 }  // end( CutPlane::ReadMult )

/*--------------------------------------------------------------------------*/

HpNum CutPlane::ReadLBMult( cIndex wFi )
{
 //??

 assert( wFi == Inf< Index >() );
 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void CutPlane::AddVariables( Index NNwVrs , cLMRow IVs )
{
 throw( NDOException( "CutPlane::AddVariables: not implemented yet." ));

 CUTLOG( 3 , std::endl << "AddVariables => "<< NNwVrs << std::endl );

 if( NumVar >= MaxNumVar )          // no space for any new variable
  return;                           // return

 if( ! NNwVrs )                     // no variables to be added
  return;                           // just return

 if( NumVar + NNwVrs > MaxNumVar )  // not enough space for all the new vars
  NNwVrs = MaxNumVar - NumVar;      // put in only the first ones

 // add the variables in the Master - - - - - - - - - - - - - - - - - - - - -

 int ncols = OsiS->getNumCols() - NrFi;

 CoinPackedVector *vectors = new CoinPackedVector[ NNwVrs ];
 for( Index i = 0 ; i < NNwVrs ; i++ )
  vectors[ i ].insert( 0 , 0 );

 for( Index n = 0 ; n < Index( OsiS->getNumRows() - 1 ) ; n++ ) {
  cIndex_Set SGBse = 0;
  cIndex SGBDm = Oracle->GetGi( Gi , SGBse, FullNames[ n ] );
  if( SGBse )
   Densify( Gi , SGBse , SGBDm, ncols + NNwVrs );

  for( Index i = 0 ; i < NNwVrs ; i++ )
   vectors[ i ].insert( n + 1 , - Gi[ i + ncols ] );
  }

 for( Index i = 0 ; i < NNwVrs ; i++ )
  OsiS->addCol( vectors[ i ] , 
                Oracle->GetUC( i ) ? - OsiS->getInfinity() : 0 ,
	        Oracle->GetUB( i ) == Inf< LMNum >() ? OsiS->getInfinity() :
		                                     Oracle->GetUB( i ) ,
		0 );
 //!!!!
 // note: the last parameter (0) in the above call is wrong, it has to be
 // the o.f. = the new entry of the 0-th subgradient

 delete[] vectors;

 if( IVs ) {
  VectAssign( Lambda + BaseSize , IVs , NNwVrs );
  Index_Set tB = Sparsify( Lambda + BaseSize , Base + BaseSize , NNwVrs ,
			   NumVar );
  if( tB > Base + BaseSize ) {
   // there actually are nonzeroes in IVs[], so the value of Fi( Lambda )
   // is no longer known
   BaseSize = tB - Base;
   *tB = Inf< Index >();

   AFiLambda = Inf< HpNum >();
   VectAssign( FiLambda , HpNum( Inf< HpNum >() ) , NrFi + 1 );
   }
  }

 NumVar += NNwVrs;

 }  // end( CutPlane::AddVariables )

/*--------------------------------------------------------------------------*/

void CutPlane::RemoveVariables( cIndex_Set whch , Index hwmny )
{
 throw( NDOException( "CutPlane::RemoveVariables: not implemented yet." ) );

 CUTLOG(3 , "\nRemoveVariables => " << hwmny-1 << std::endl );

 if( ! whch ) {
  int* vars = new int[ NumVar ];
  for( Index i = 0 ; i < NumVar ; i++ )
   vars[ i ] = ( i + NrFi );
  OsiS->deleteCols( NumVar , vars );
  delete[] vars;
  }
 else {
  int* vars = new int[ hwmny - 1 ];
  for( Index i = 0 ; i < hwmny - 1 ; i++ )
   vars[ i ] = whch[ i ] + NrFi;
  OsiS->deleteCols( hwmny - 1 , vars );
  delete[] vars;
  }

 // Delete removed Lambdas
 if( Base )
  Densify( Lambda , Base , BaseSize , NumVar);

 if( whch ) {
  for( Index i = 0 ; i < hwmny - 1 ; i++ )
   Lambda[ whch[ i ] ] = 0;
  }
 else
  for( Index i = 0 ; i < NumVar ; i++ )
   Lambda[ i ] = 0;
  
 NumVar = Oracle->GetNumVar();
 if( Base )
  delete[] Base;

 Base = new Index[ NumVar ];
 Index_Set BaseEnd;
 BaseEnd = Sparsify( Lambda , Base , NumVar );

 if( Base[ NumVar ] != BaseEnd[ 0 ] ) {  // Lambda is sparse
  BaseEnd[ 0 ] = Inf< Index >();
  for( BaseSize = 0 ; Base[ BaseSize ] != Inf< Index >() ; BaseSize++ );
  }
 else {  // Lambda is dense
  delete[] Base;
  Base = 0;
  BaseSize = NumVar;
  }

 Oracle->SetLamBase( Base , BaseSize );
 Oracle->SetLambda( Lambda );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

CutPlane::~CutPlane()
{
 // memory deallocation - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Oracle )
  MemDealloc();

 }  // end( ~CutPlane )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

inline void CutPlane::FiAndGi( void )
{
 // pass Lambda to the Oracle - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Oracle->SetLamBase( BaseSize < NumVar ? Base : 0 , BaseSize );
 Oracle->SetLambda( Lambda );

 // compute all Fi( i ) values- - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index FiSp = 1 ; FiSp <= NrFi ; FiSp++ ) {
  FiLambda[ FiSp ] = Oracle->Fi( FiSp );
  FiEvaltns++;

  if( FiLambda[ FiSp ] == Inf< HpNum >() ) {
   CUTLOG( 1 , "\tError: Fi( Lambda ) = - INF not supported." << std::endl );
   Result = kError;
   return;
   }

  if( FiLambda[ FiSp ] == - Inf< HpNum >() ) {
   CUTLOG( 1 , "\tFi( Lambda ) = INF => STOP." << std::endl );
   Result = kUnbndd;
   return;
   }
  }  // end( for( FiSp ) )

 // compute AFiLambda, update ABestFiLambda - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiLambda[ 0 ] = Oracle->Fi( 0 );
 AFiLambda = SumV( FiLambda , NrFi + 1 );

 if( AFiLambda < ABestFiLambda ) {
  ABestFiLambda = AFiLambda;
  VectAssign( BestFiLambda , FiLambda , NrFi + 1 );

  if( KpBstL ) {
   if( BaseSize < NumVar ) {
    cLMRow tL = Lambda;
    cIndex_Set tB = Base;
    for( Index i = 0 ; i < NumVar ; i++ )
     if( *tB == i ) {
      BestLambda[ i ] = *(tL++);
      tB++;
      }
     else
      BestLambda[ i ] = 0;
    }
   else
    VectAssign( BestLambda , Lambda , NumVar );
   }
  }

 // now start collecting subgradients - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // setting parameters for the heuristic cycle

 Index HeurStop , HeurInt;
 Index HeurIter = 0;
 if( AValueRMP <= LowerBound ) {  // when the LB is tight, it means that
  HeurStop = HeurSubGU;           // the RMP is "artificially" bounded
  HeurInt = 1;
  }
 else {                           // the RMP is "truly" bounded
  HeurStop = HeurSubGF;
  HeurInt = HeurInvl;
  }

 do {
  // begin subproblems cycle - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index FiSp = 1 ; FiSp <= NrFi ; FiSp++ ) {
   // test if the row is effective - - - - - - - - - - - - - - - - - - - - -

   CUTLOG( 3 , std::endl << "Gi[ " << FiSp << " ]: " );

   if( ( ValueRMP[ FiSp ] < - Inf< HpNum >() ) &&
       ( FiLambda[ FiSp ] - ValueRMP[ FiSp ] <=
	 EpsLin * std::max( ABS( AFiLambda ) , HpNum( 1 ) ) / NrFi ) ) {
     CUTLOG( 3 , "not effective." );
    continue;  // it is not: loop
    }

   if( ! Oracle->NewGi( FiSp ) ) {  // there is no new subgradient
    CUTLOG( 3 , "not available." );
    continue;  // loop
    }

   RowsAdded++;  // something new and useful has been found
   CUTLOG( 3 , "added." );

   // fetch the subgradient from the Oracle- - - - - - - - - - - - - - - - -

   cIndex_Set SGBse;
   cIndex SGBDm = Oracle->GetGi( Gi , SGBse );

   GiEvaltns++;
   HpNum eps = Oracle->GetVal();

   // expand array of names if it is full- - - - - - - - - - - - - - - - - -

   if( KpPrimals ) {
    if( FirstEmpty == NumNames * AllocNumber ) {
     FullNames = (Index*) realloc( FullNames , 
				   sizeof( Index ) * NumNames *
				   ( AllocNumber + 1 ) );
     if( ! FullNames )
      throw( NDOException( "Error in memory reallocation" ) );

     SubPNames = (Index*) realloc( SubPNames ,
				   sizeof( Index ) * NumNames *
				   ( AllocNumber + 1 ) );
     if( ! SubPNames )
      throw( NDOException( "Error in memory reallocation" ) );

     EmptyNames = (Index*) realloc( EmptyNames ,
				    sizeof(Index) * NumNames *
				    ( AllocNumber + 1 ) );
     if( ! EmptyNames )
      throw( NDOException( "Error in memory reallocation" ) );

     for( Index i = NumNames * AllocNumber ;
	  i < NumNames * ( AllocNumber + 1 ) ; i++ )
      EmptyNames[ i ] = i;

     AllocNumber++;
     Oracle->SetMaxName( NumNames * AllocNumber );
     }

    // tell the name of the subgradient to the FiOracle - - - - - - - - - - -

    Oracle->SetGiName( EmptyNames[ FirstEmpty ] );
    FullNames[ LastFull + 1 ] = EmptyNames[ FirstEmpty ];
    FirstEmpty++;
    LastFull++;
    SubPNames[ LastFull ] = FiSp;

    }  // end( if( KpPrimals ) )

   // add row to the RMP - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // Gi is an eps-subgradient of the FiSp component of Fi() at Lambda, which
   // means that the inequality
   //   v[ FiSp ] >= Fi( FiSp ) + Gi * ( L - Lambda ) - eps
   // is valid for each L. Hence, the constraint to be added to the RMP is
   //   v[ FiSp ] - Gi * L >= Fi( FiSp ) - Gi * Lambda - eps

   double rhs = double( FiLambda[ FiSp ] ) - eps;

   CoinPackedVector vect;
   vect.insert( FiSp - 1 , 1 );

   if( SGBse ) {             // Gi is sparse
    for( Index i = 0 ; i < SGBDm ; i++ )
     vect.insert( SGBse[ i ] + NrFi , - Gi[ i ] );

    if( BaseSize < NumVar )  // Lambda is also sparse
     rhs -= ScalarProduct( Gi , SGBse , Lambda , Base );
    else                     // Lambda is dense
     rhs -= ScalarProduct( Gi , Lambda , SGBse );
    }
   else {                    // Gi is dense
    for( Index i = 0 ; i < NumVar ; i++ )
     vect.insert( i + NrFi , - Gi[ i ] );

    if( BaseSize < NumVar )  // but Lambda is sparse
     rhs -= ScalarProduct( Lambda , Gi , Base );
    else                     // both are dense
     rhs -= ScalarProduct( Lambda , Gi , NumVar );
    }

   OsiS->addRow( vect , 'G' , rhs , 0 );

   }  // end( for( FiSp ) )
  } while( ( ++HeurIter < HeurStop ) && ( ParIter % HeurInt == 0 ) );

 // try to update the LowerBound- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cHpNum LwrBnd = Oracle->GetLowerBound();
 if( LwrBnd > LowerBound ) {
  LowerBound = LwrBnd;
  TrueLB = true;
  OsiS->setRowLower( 0 , LwrBnd );
  }
 }  // end( CutPlane::FiAndGi )

/*--------------------------------------------------------------------------*/

inline Index CutPlane::DelRMPRows( HpNum lim )
{
 cIndex count = OsiS->getNumRows();

 const double *rhsCol = OsiS->getRightHandSide();
 const double *axCol = OsiS->getRowActivity();
 HpRow slacks = new HpNum[ count ];

 for( Index i = 0 ; i < count ; i++ )
  slacks[ i ] = rhsCol[ i ] - axCol[ i ];

 Index delCount = 0 , newCount = 0;
 int *delRows = new int[ count + 1 ];

 if( KpPrimals ) {  // keep track of names- - - - - - - - - - - - - - - - - -
  Index_Set newNames = new Index[ NumNames * AllocNumber ];
  Index_Set newSPNames = new Index[ NumNames * AllocNumber ];
  for( Index i = 1 ; i < count ; i++ ) {
   if( slacks[ i ] < -lim )
    delRows[ delCount++ ] = i;
   else {
    newSPNames[ newCount ] = SubPNames[ i - 1 ];
    newNames[ newCount++ ] = FullNames[ i - 1 ];
    // [ i - 1 ] because the first row of the RMP does not havea name
    }
   }

  if( delCount ) {
   OsiS->deleteRows( delCount , delRows );
   for( Index i = 0 ; i < delCount ; i++ ) {
    FirstEmpty--;
    LastFull--;
    EmptyNames[ FirstEmpty ] = FullNames[ delRows[ i ] - 1 ];
    Oracle->Deleted( FullNames[ delRows[ i ] - 1 ] );
    }

   delete[] FullNames;
   FullNames = newNames;
   delete[] SubPNames;
   SubPNames = newSPNames;
   }
  else
   delete[] newNames;
  }
 else {  // don't bother with names - - - - - - - - - - - - - - - - - - - - -
  for( Index i = 1 ; i < count ; i++ )
   if( slacks[ i ] < - lim )
    delRows[ delCount++ ] = i;

  if( delCount )
   OsiS->deleteRows( delCount , delRows );
  }
  
 delete[] delRows;
 delete[] slacks;

 return( delCount );

 }  // end( CutPlane::DelRMPRows )

/*--------------------------------------------------------------------------*/

inline void CutPlane::OutRsts( void )
{
 // output results, times and statistics - - - - - - - - - - - - - - - - - - -

 CUTLOG( 1 , std::endl << "Solve: " << SCalls );
 CUTLOG( 1 , std::endl << "Solution Value: " << AValueRMP );
 CUTLOG2( 1 , NDOt , std::endl << "Tot. time (sec): " << NDOt->Read() );
 CUTLOG( 1 , std::endl << "Gap: " << Gap );
 CUTLOG( 1 , std::endl << "Iterations: " << ParIter );
 CUTLOG( 1 , std::endl << "Fi() time (sec): " << Oracle->FiTime() );
 CUTLOG( 1 , std::endl << "Max Number of rows (RMP): " << MaxNRows );
 CUTLOG( 1 , std::endl << "Final Number of rows (RMP): "
	               << OsiS->getNumRows() );
 CUTLOG( 1 , std::endl << "Final Number of columns (RMP): "
	               << OsiS->getNumCols() );
 CUTLOG( 1 , std::endl << "Total Fi() evaluations : " << FiEvaltns
	               << std::endl );
 }

/*--------------------------------------------------------------------------*/

inline void CutPlane::MemDealloc( void )
{
 delete[] Gi;
 delete[] Lambda;
 delete[] BestLambda;
 delete[] BestFiLambda;
 delete[] FiLambda;
 delete[] ValueRMP;
 delete[] FullNames;
 delete[] EmptyNames;
 delete[] SubPNames;

 }  // end( MemDealloc )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File CutPlane.C  --------------------------*/
/*--------------------------------------------------------------------------*/
