/*--------------------------------------------------------------------------*/
/*----------------------------- File Bundle.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the Bundle class, which implements the NDOSolver
 * interface for NonDifferentiable Optimization Solvers, as described
 * NDOSlver.h, using a "Generalized Bundle" algorithm.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>

#include "OPTvect.h"

#include "Bundle.h"

#include <assert.h>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define NOISE_REDUCTION_FIRST 1

/* If NOISE_REDUCTION_FIRST > 0, then in case of trouble the Bundle will first
   try to raise t, then ask the FiOracle to increase the accuracy. Otherwise
   the order is reversed. */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#define CHECK_DS 0

/* If CHECK_DS > 0, various data structures are checked for correctness
   during the run of the algorithm, tremendously slowing down the algorithm
   but allowing to debug the thing.

   What data structures are checked is coded bit-wise in CHECK_DS:

    bit 0 (+ 1)  =>  the value of linearization errors of subgradients /
                     RHS of constraints as kept updated within the Master
		     Problem is compared with that computed by the FiOracle
		     [see GetVal() in FiOracle.h]; this may be useful while
		     debugging a new MPSolver or a new FiOracle.

   CHECK_DS > 0 forces asserts() to work within this unit. */

#if( CHECK_DS )
 #ifdef NDEBUG
  #undef NDEBUG
 #endif
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( LOG_BND )
 #define BLOG( l , x ) if( NDOLLvl > l ) *NDOLog << x

 #define BLOG2( l , c , x ) if( ( NDOLLvl > l ) && c ) *NDOLog << x

 #define BLOGb( l , x ) if( NDOLLvl & l ) *NDOLog << x

 #define BLOG2b( l , c , x ) if( ( NDOLLvl & l ) && c ) *NDOLog << x
#else
 #define BLOG( l , x )

 #define BLOG2( l , c , x )

 #define BLOGb( l , x )

 #define BLOG2b( l , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static const HpNum Nearly  = 1.01;
static const HpNum Nearly2 = 1.02;

static const HpNum DefMPEFsb = 1e-6;  // default value for MPEFsb
static const HpNum DefMPEOpt = 1e-6;  // default value for MPEOpt

static const char LogBnd = 16;        // log Bundle changes
static const char LogVar = 32;        // log variables changes

static cIndex tSP1Msk = ~ 3;          // mask for tSPar1
static cIndex kSLTTS =  4;            // "soft" long-term t-strategy
static cIndex kHLTTS =  8;            // "hard" long-term t-strategy
static cIndex kBLTTS = 12;            // "balancing" long-term t-strategy
static cIndex kEGTTS = 16;            // "endgame" long-term t-strategy

static const unsigned char RstAlg =  1;  // don't reset algorithmic parameters
static const unsigned char RstCrr =  2;  // don't reset current point
static const unsigned char RstSbg =  4;  // don't reset subgradients
static const unsigned char RstCnt =  8;  // don't reset constraints
static const unsigned char RstFiV = 16;  // don't reset FiVals

static cIndex InINF = Inf< Index >();
static cHpNum HpINF = Inf< HpNum >();

/*--------------------------------------------------------------------------*/
/*---------------------- IMPLEMENTATION OF Bundle --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

Bundle::Bundle( std::istream *iStrm )
        :
        NDOSolver( iStrm )
{
 Master = NULL;

 // initialize algorithmic parameters - - - - - - - - - - - - - - - - - - - -

 DfltdSfInpt( iStrm , BPar1 , int( 10 ) );
 DfltdSfInpt( iStrm , BPar2 , int( 100 ) );
 DfltdSfInpt( iStrm , BPar3 , HpNum( -1 ) );
 DfltdSfInpt( iStrm , BPar4 , HpNum( -1 ) );
 DfltdSfInpt( iStrm , BPar5 , HpNum( 30 ) );
 DfltdSfInpt( iStrm , BPar6 , int( 0 ) );

 DfltdSfInpt( iStrm , m1 , HpNum( .1 ) );
 DfltdSfInpt( iStrm , m3 , HpNum( 3 ) );

 DfltdSfInpt( iStrm , mxIncr , HpNum( 10 ) );
 DfltdSfInpt( iStrm , mnIncr , HpNum( 1.5 ) );
 DfltdSfInpt( iStrm , MnSSC , int( 0 ) );
 DfltdSfInpt( iStrm , mxDecr , HpNum( .1 ) );
 DfltdSfInpt( iStrm , mnDecr , HpNum( .66 ) );
 DfltdSfInpt( iStrm , MnNSC , int( 0 ) );

 if( mxIncr < 1 ) mxIncr = 1;
 if( mnIncr < 1 ) mnIncr = 1;
 if( mxDecr > 1 ) mxDecr = 1;
 if( mnDecr > 1 ) mnDecr = 1;

 DfltdSfInpt( iStrm , tMaior , HpNum( 1e+6 ) );
 DfltdSfInpt( iStrm , tMinor , HpNum( 1e-6 ) );
 DfltdSfInpt( iStrm , tInit  , HpNum( 1 ) );

 DfltdSfInpt( iStrm , tSPar1 , int( 0 ) );
 DfltdSfInpt( iStrm , tSPar2 , HpNum( .1 ) );

 DfltdSfInpt( iStrm , PPar1 , int( 30 ) );
 DfltdSfInpt( iStrm , PPar2 , int( 10 ) );
 DfltdSfInpt( iStrm , PPar3 , int( 5 ) );

 DfltdSfInpt( iStrm , MPEFsb , DefMPEFsb );
 DfltdSfInpt( iStrm , MPEOpt , DefMPEOpt );

 // some initializations- - - - - - - - - - - - - - - - - - - - - - - - - - -

 t = tInit;
 Prevt = HpINF;
 KpBstL = false;
 tHasChgd = LHasChgd = true;

 IsEasy = 0;
 Lambda = 0;
 Lambda1 = 0;
 LmbdBst = 0;
 LamBase = 0;
 Lam1Bse = 0;
 FiLambda = 0;
 FiBest = 0;
 FiLambda1 = 0;
 RfrncFi = 0;
 whisZ = 0;
 whisG1 = 0;
 ScPr1 = 0;
 Alfa1 = 0;
 DeltaAlfa = 0;
 LowerBound = 0;
 FreList = 0;
 OOBase = 0;
 InctvCtr = 0;
 nBase = 0;

 #if( NONMONOTONE )
  FiVals = 0;
 #endif
 
 }  // end( Bundle )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void Bundle::SetMPSolver( MPSolver *MPS )
{
 Master = MPS;
 if( Master && Oracle )
  InitMP();

 }  // end( Bundle::SetMPSolver )

/*--------------------------------------------------------------------------*/

void Bundle::SetFiOracle( FiOracle *Fi )
{
 if( Oracle ) {  // changing from a previous oracle - - - - - - - - - - - - -
                 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  RemoveItems();  // clear the bundle
  MemDealloc();   // deallocate memory
  Oracle = NULL;  // discard all references to old oracle
  }

 if( Fi ) {  // setting a new oracle- - - - - - - - - - - - - - - - - - - - -
             // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // "throw" the method of the base class - - - - - - - - - - - - - - - - - -

  NDOSolver::SetFiOracle( Fi );

  // read information about Fi- - - - - - - - - - - - - - - - - - - - - - - -

  MaxNumVar = Oracle->GetMaxNumVar();

  if( NrFi > 1 ) {
   IsEasy = new bool[ NrFi ];
   IsEasy--;
   bool HasEasy = false;

   for( Index k = 0 ; k++ < NrFi ; )
    if( Oracle->GetBNC( k ) )
     IsEasy[ k ] = HasEasy = true;
    else
     IsEasy[ k ] = false;

   if( ! HasEasy ) {
    delete[] ++IsEasy;
    IsEasy = NULL;
    }
   }
  else
   IsEasy = NULL;

  // allocate memory- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( PPar2 ) {  // L.V.G. "on"
   LamBase = new Index[ MaxNumVar + 3 ];  // the set of "active" variables
   nBase   = new Index[ MaxNumVar + 3 ];  // temporay for changes in LamBase

   *(LamBase++) = InINF;  // set a "floor" to the bases: this is used
   *(nBase++) = InINF;    // in FormD() when they are read backwards

   LamBase[ MaxNumVar + 1 ] = InINF;  // put a stop to the second part
   nBase[ MaxNumVar + 1 ] = InINF;    // of the bases, which will hold
                                             // the indices of the variables
                                             // added/removed

   LamBase[ LamDim = 0 ] = InINF;  // at start no variables are active

   Lam1Bse = new Index[ MaxNumVar + 1 ];  // "active" variables for Lambda1
   *Lam1Bse = InINF;

   if( PPar3 )                            // variables are also deleted
    InctvCtr = new Index[ MaxNumVar ];    // counter for removing variables
   else {
    InctvCtr = NULL;
    PPar1 = InINF;
    }
   }
  else {         // L.V.G. "off"
   LamBase = nBase = Lam1Bse = NULL;
   LamDim = NumVar;
   InctvCtr = NULL;
   }

  Lambda = new LMNum[ MaxNumVar ];             // the current point
  VectAssign( Lambda , LMNum( 0 ) , NumVar );  // the default starting point

  Lambda1 = new LMNum[ MaxNumVar ];            // the tentative point
  if( KpBstL )
   LmbdBst = new LMNum[ MaxNumVar ];           // best point found so far
  else
   LmbdBst = NULL;

  OOBase = new SIndex[ BPar2 ];
  VectAssign( OOBase , SIndex( Inf< SIndex >() ) , BPar2 );
  // counter for eliminating outdated items: Inf< SIndex >() means empty
  FreList = new Index[ BPar2 ];                // list of free bundle slots
  whisZ   = new Index[ NrFi ];                 // for each component, the
  whisZ--;                                     // name of its "Z" if it is
                                               // in the bundle
  FiLambda  = new HpNum[ NrFi + 1 ];           // current, ...
  FiLambda1 = new HpNum[ NrFi + 1 ];           // tentative, ...
  FiBest    = new HpNum[ NrFi + 1 ];           // best, ...
  RfrncFi   = new HpNum[ NrFi + 1 ];           // and reference Fi() values

  FiStatus = new FiOracle::FiStatus[ NrFi ];
  VectAssign( FiStatus , FiOracle::kFiStop , NrFi );
  FiStatus--;

  LowerBound = new HpNum[ NrFi + 1 ];          // lower bounds
  VectAssign( LowerBound , - HpINF , NrFi + 1 );
  TrueLB = false;

  VectAssign( RfrncFi , HpNum( 0 ) , NrFi + 1 );
  *FiLambda1 = *FiBest = HpINF;         // Fi( Lambda ) is not known

  whisG1 = new Index[ NrFi ];
  VectAssign( whisG1 , Index( InINF ), NrFi );
  // no representative yet
  whisG1--;

  ScPr1 = new LMNum[ NrFi + 1 ];
  Alfa1 = new HpNum[ NrFi + 1 ];
  VectAssign( ScPr1 , HpNum( 0 ) , NrFi + 1 );
  VectAssign( Alfa1 , HpNum( 0 ) , NrFi + 1 );
  DeltaAlfa = new HpNum[ NrFi ];
  DeltaAlfa--;

  FreDim = 0;
  BHasChgd = true;  // ensure SetLamBase() is called at least once
  Result = kError;
  SSDone = false;

  ReSetAlg( RstCrr | RstSbg | RstCnt );       // Fi( Lambda ) is reset inside

  // tell the oracle about the NDOSolver and its settings - - - - - - - - - -

  Oracle->SetNDOSolver( this );
  Oracle->SetMaxName( BPar2 );
  Oracle->SetPrecision( EpsFi = ABS( EInit  ) );

  // warning: the following things can only be done *after* that
  // Oracle->SetMaxName() has been invoked, because they use methods of the
  // oracle which depends on knowledge of the MaxName to work properly
  // read b0- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // here one could initialize b0, if that was found to be of any use
  // b0 = Oracle->GetVal( BPar2 );

  // initialize the MP Solver, if any - - - - - - - - - - - - - - - - - - - -

  if( Master )
   InitMP();
  }
 else  // Fi == NULL- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( Master )        // a MP solver is set
   Master->SetDim();  // clear all its internal state

 }  // end( Bundle::SetFiOracle )

/*--------------------------------------------------------------------------*/

void Bundle::SetLambda( cLMRow tLambda )
{
 if( ! Oracle )
  throw( NDOException( "Bundle::SetLambda: Oracle == NULL" ) );

 // set Lambda1 = NewLambda - Lambda- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( tLambda ) {  // ... if a Lambda is given - - - - - - - - - - - - - - - -
  VectAssign( Lambda1 , tLambda , NumVar );

  if( PPar2 )
   VectSubtractB( Lambda1 , Lambda , LamBase );
  else
   VectSubtract( Lambda1 , Lambda , NumVar );
  }
 else           // tLambda = 0- - - - - - - - - - - - - - - - - - - - - - - -
  if( PPar2 )
   VectMAssignB( Lambda1 , Lambda , LamBase , NumVar );
  else
   VectMAssign( Lambda1 , Lambda , NumVar );

 // set Lambda- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( PPar2 ) {  // sparse Lambda- - - - - - - - - - - - - - - - - - - - - - -

  LamDim = 0;
  if( tLambda ) {
   cLMNum teD = Master->EpsilonD();
   for( Index i = 0 ; i < NumVar ; i++ ) {
    cLMNum Li = *(tLambda++);
    if( ABS( Li ) > teD ) {
     Lambda[ LamDim ] = Li;
     LamBase[ LamDim++ ] = i;
     if( PPar3 )
      InctvCtr[ i ] = 0;
     }
    }
   }

  LamBase[ LamDim ] = InINF;
  BHasChgd = true;

  Master->SetActvSt( LamBase , LamDim );  // re-set the active set
  }
 else {       // dense Lambda - - - - - - - - - - - - - - - - - - - - - - - -
  if( tLambda )
   VectAssign( Lambda , tLambda , NumVar );
  else
   VectAssign( Lambda , LMNum( 0 ) , NumVar );

  LamDim = NumVar;  // keep LamDim == NumVar
  }

 // change the Current Point in the subproblem solver - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // pretend that the value of Fi() in the *new* Lambda is *the same as in the
 // old Lambda*, as the correction will be done as soon as Fi() is computed

 VectAssign( FiLambda1 , HpNum( 0 ) , NrFi + 1 );
 *FiLambda = *FiBest = HpINF;  // Fi( Lambda ) is not known

 Master->ChangeCurrPoint( Lambda1 , FiLambda1 );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SSDone = false;

 }  // end( SetLambda )

/*--------------------------------------------------------------------------*/

void Bundle::KeepBestLambda( const bool KBL )
{
 if( KpBstL != KBL ) {
  if( KpBstL ) {
   delete[] LmbdBst;
   LmbdBst = NULL;
   }
  else
   LmbdBst = new LMNum[ MaxNumVar ];

  KpBstL = KBL;
  }
 }  // end( Bundle::KeepBestLambda )

/*--------------------------------------------------------------------------*/

void Bundle::SetNDOLog( std::ostream *outs , const char lvl )
{
 NDOSolver::SetNDOLog( outs , lvl );

 #if( LOG_BND )
  if( NDOLLvl > 1 )
   *NDOLog << std::endl << "Vars = " << NumVar << " (" << MaxNumVar
	   << ") ~ Max # = " << BPar2 << " ~ Rfrsh = " << BPar1 << " ~ t* = "
	   << tStar << " ~ Eps = " << EpsLin << std::endl
	   << "t = " << tInit << " in [" << tMinor << ", " << tMaior
           << "] ~ Incr = " << mxIncr << " | " << mnIncr << " | " << MnSSC
	   << " ~ Decr = " << mxDecr << " | " << mnDecr << " | " << MnNSC
	   << std::endl
	   << "m1 = " << m1 << " ~ m3 = " << m3 << " ~ tS = " << tSPar1
           << "(" << tSPar2 << ") ~ Pricing: " << PPar1 << " - " << PPar2
	   << " - " << PPar3 << std::endl
	   << "FiEps = " << EpsCurr  << " in [ " << EInit << " , "
	   << EFnal << " ] with decr = " << EDcrs << " (" << EStps << ")"
	   << std::endl;
 #endif

 }  // end( Bundle::SetNDOLog )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

NDOSolver::NDOStatus Bundle::Solve( void )
{
 // basic sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! Master )
  throw( NDOException( "Bundle::Solve: Master not set yet" ) );

 if( ! Oracle )
  throw( NDOException( "Bundle::Solve: Oracle not set yet" ) );

 // initializations - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NDOt )
  NDOt->Start();

 Result = kOK;
 SCalls++;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle starts here- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 do {
  // construct the direction d- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormD();

  if( Result )  // problems in the Master Problem solver
   break;

  // a little bookkeeping - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // update out-of-base counters- - - - - - - - - - - - - - - - - - - - - - -

  UpdtCntrs();

  // some log - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Log1();

  // contrast MPsolver::Alfa[] with Oracle::Alfa[]- - - - - - - - - - - - - -

  StrongCheckAlfa();

  // hook for derived classes - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  EIStatus rs = EveryIteration();

  if( rs == kEIAbort ) {
   BLOG( 1 , " EveryIteration():STOP" << std::endl );
   Result = kStopped;
   break;
   }

  if( rs == kEILoopNow ) {
   BLOG( 1 , " EveryIteration():loop" << std::endl );
   continue;
   }

  // check the status of the FiOracle and take the necessary action - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FiOracle::FiStatus fs = Oracle->GetFiStatus();

  if( fs == FiOracle::kFiStop ) {
   BLOG( 1 , " ~ FiOracle:STOP" << std::endl );
   Result = kStopped;
   break;
   }

  if( fs == FiOracle::kFiError ) {
   BLOG( 1 , " ~ Error in the FiOracle" << std::endl );
   Result = kError;
   break;
   }

  if( fs == FiOracle::kFiChgd ) {
   BLOG( 1 , " ~ Fi changed: loop" << std::endl );
   continue;
   }

  // check for optimality - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( rs != kEIContAnyway ) && ( fs != FiOracle::kFiCont ) )
   if( IsOptimal() )
    break;

  // check for running time - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( MaxTime && NDOt )
   if( NDOt->Read() > MaxTime ) {
    Result = kStpTime;
    break;
    }

  // Hard Long-Term t-strategy- - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the hard long-term t-strategy requires t to increase if the step is too
  // small, and therefore has to be checked before the others

  if( ( ( tSPar1 & tSP1Msk ) == kHLTTS ) && ( *FiLambda < HpINF ) ) {
   HpNum AFL = ABS( *FiLambda );
   if( AFL < 1 )
    AFL = 1;

   if( vStar <= tSPar2 * EpsU * AFL ) {
    BLOG( 1 , "small v => increase t" << std::endl << "           " );

    // collect two numbers vc and vl such that v( tNew ) >= vc + tNew * vl
    // we require that v( tNew ) >= vc + tNew * vl = tSPar2 * EpsU * AFL
    // ==> tNew = ( tSPar2 * EpsU * AFL - vc ) / vl

    HpNum vl , vc;
    Master->SensitAnals( vl , vc );

    HpNum tt;
    if( - vl < Eps< HpNum >() )  // v( t ) is [almost] constant ==> D*_t [~]= 0
     tt = tStar;       // ==> the CP model is ~bounded
    else
     tt = std::min( tStar , ( tSPar2 * EpsU * AFL * Nearly + vc ) /
		            ( - vl ) );

    if( ( tHasChgd = ( tt != t ) ) ) {
     t = tt;
     continue;         // loop only if t changes
     }
    }
   }  // end if( Hard t-strategy )

  // a real iteration (iterations where Fi() is not evaluated do not count) -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ParIter++;

  // change the "precision" in computing Fi() - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is the "regular mechanism"; gradually decrease the precision along
  // the iterations

  if( EStps ) {
   cIndex nstps = ::ceil( double( EStps > 0 ? NrIter() : NrSSs() ) /
			  double( ABS( EStps ) ) );
   EpsCurr = EDcrs >= 0 ? ABS( EInit ) : EpsU;
   if( EFnal >= 0 )
    EpsCurr *= pow( ABS( EDcrs ) , EFnal * nstps );
   else
    EpsCurr *= ABS( EDcrs ) * ( nstps ? pow( nstps , EFnal ) : 1 );

   EpsCurr = std::max( EpsLin , std::min( EpsCurr , ABS( EInit ) ) );
   }
  else
   EpsCurr = EDcrs >= 0 ? ABS( EInit ) : ( ABS( EDcrs ) * EpsU );

  if( ( EpsCurr < EpsFi ) ||
      ( ( EpsCurr > EpsFi ) && ( EInit < 0 ) ) ) {
   // only allow increasing EpsFi if EInit < 0, always allow decreasing it
   BLOG( 1 , " ~ changing precision to " << EpsCurr << std::endl
	 << "           " );
   Oracle->SetPrecision( EpsFi = EpsCurr );
   }

  // calculate Lambda1- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  FormLambda1( t );

  // update the number of items to be fetched from the oracle - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  UpdtaBP3();

  // eliminate outdated info- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // This is done *after* the call to Oracle->GetFiStatus(), as well as after
  // the call to Master->SensitAnals() in the Hard Long-Term t-strategy and
  // to FormLambda1(), because elimination of items from the bundle may make
  // the current solution of the master problem invalid, and therefore all
  // solution information may be lost. In theory this should not happen, since
  // only items "out of base" are eliminated, and therefore the solution
  // remains optimal; however, not all MPSolvers may behave in this respect.

  SimpleBStrat();

  // calculate Fi( Lambda1 )- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  bool MPchgs;  // true if no cycling will occur
  for( ;; ) {   // ... possibly more than once due to precision issues

   // actually compute Fi and collect subgrads- - - - - - - - - - - - - - - -
   // meanwhile check if the new subgradients change the CP model enough
   MPchgs = FiAndGi();

   if( Result == kError )
    break;

   // if there are negative alphas, something has indeed changed
   if( *FiLambda < HpINF )
    MPchgs |= CheckAlfa();

   MPchgs |= DoSS();  // doing a SS clearly changes the MP

   if( MPchgs )       // if something changes
    break;            // all done

   // check for running time - - - - - - - - - - - - - - - - - - - - - - - - -
   // if we get here there is something wrong with the FiOracle's precision;
   // the possible solution is to give it a few more resources to try to do
   // it, but this is not possible if we have ran out of time

   if( MaxTime && NDOt )
    if( NDOt->Read() > MaxTime ) {
     Result = kStpTime;
     break;
     }

   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // if neither a SS is done nor the Master Problem changes, there must be
   // something wrong with the FiOracle's precision

   #if( NOISE_REDUCTION_FIRST  )
    // we try to patch this by playing with t first, which corresponds to
    // saying that we think it faster to re-solve the Master Problem rather
    // than to get more precision from the FiOracle

    if( ( DSTS >= EpsLin * std::max( HpNum( 1 ) , ABS( *FiLambda ) ) ) &&
	( t < tMaior ) ) {
     t = std::min( t * mxIncr , tMaior );
     BLOG( 1 , " ~ noise reduction: t increased to " << t << std::endl );
     tHasChgd = true;
     break;            
     }
   #endif

   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // Now check if any of the components has not been solved with the required
   // accuracy yet; in this case allow the component(s) to be re-computed
   // again (to increase accuracy). Note that you expect this to happen only
   // a limited number of times, both because eventually a "good enough"
   // solution will be found by the oracle, and because eventually time will
   // run out.

   FiOracle::FiStatus FiStt = FiOracle::kFiNorm;
   for( Index k = 0 ; k++ < NrFi ; )
    if( ( ( ! IsEasy ) || ( ! IsEasy[ k ] ) ) &&
	( FiStatus[ k ] == FiOracle::kFiStop ) ) {
     FiStt = FiOracle::kFiStop; 
     break;
     }

   if( FiStt == FiOracle::kFiStop )  // at least one non-easy component did
    continue;                        // not finish, go give them more time

   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // the oracle reports to have computed the function to the required
   // precision, but this is not enough: here goes the "emergency mechanism"
   // that tries to increase the accuracy, but of course this all depends on
   // if the oracle is actually available to do it

   EpsCurr = std::max( EpsLin , std::min( EpsCurr , EpsFi ) * ABS( EDcrs ) );

   if( EpsCurr < EpsFi )  // unless we have already hit EpsLin
    if( Oracle->SetPrecision( EpsFi = EpsCurr ) ) {  // if it is
     BLOG( 1 , " ~ increasing precision to " << EpsCurr << std::endl );
     VectAssign( FiStatus + 1 , FiOracle::kFiStop , NrFi );
     // reset FiStatus[]: if a component has said "kFiOK" before, with a
     // coarser accuracy, this does not mean the computation is still OK 
     // now that the accuracy has increased, so we must assume it is not
     continue;                                       // go compute Fi() again
     }

   #if( ! NOISE_REDUCTION_FIRST  )
    // we try to patch this by playing with t last, which corresponds to
    // saying that we think it faster to get more precision from the
    // FiOracle rather than to re-solve the Master Problem 

    if( ( DSTS >= EpsLin * std::max( HpNum( 1 ) , ABS( *FiLambda ) ) ) &&
	( t < tMaior ) ) {
     t = std::min( t * mxIncr , tMaior );
     BLOG( 1 , " ~ noise reduction: t increased to " << t << std::endl );
     tHasChgd = true;
     break;            
     }
   #endif

   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // Now the Bundle is potentially in trouble: the new informations does not
   // change the model, changing t does not help, and the oracle is not
   // willing to provide any more precision right now. However, the oracle
   // may "know" this, this being an iteration where the Bundle would have
   // stopped already, but it was told not to by the oracle. This is clearly
   // at risk of cycling, so the oracle must have some internal mechanism to
   // avoid that. in this case, trust the oracle: just go to the next
   // iteration and hope that something eventually will change.

   if( fs == FiOracle::kFiCont ) {  // if the oracle wants to rule
    MPchgs = true;                  // pretend to believe no cycling will
    break;                          // occur since the oracle says so
    }
   else {                           // really, nothing else to do but quit
    BLOG( 1 , " ~ too low precision in the FiOracle" << std::endl );
    Result = kLwPrcsn;
    break;
    }
   }  // end( for( ever ) )

  // check whether the Lower Bounds have changed- - - - - - - - - - - - - - -

  UpdtLowerBound();

  // some log about the newly obtained information- - - - - - - - - - - - - -

  Log2();

  if( *FiLambda1 == - HpINF ) {
   Result = kUnbndd;
   break;
   }

  // check whether either any error has occurred or time has expired- - - - -

  if( ( Result == kError ) || ( Result == kStpTime ) )
   break;

  if( Result == kLwPrcsn )
   break;

  if( ( ! MPchgs ) && tHasChgd )  // "noise reduction": t has changed,
   continue;                      // so go solve the master problem again
                                  // (no NS/SS decision can be made)

  // check the Lower Bound- - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // note: TrueLB is true if *LowerBound > GetMinusInfinity(). Termination
  //       "by objective function value only" is only enabled if TrueLB is
  // true. In the Lagrangian case, if one sets as LowerBound the value of a
  // feasible solution it may stop here is that solution is EpsLin-optimal.
  // However, doing so might "disrupt the convexified solution", because the
  // Master Problem is not solved and therefore the optimal multipliers are
  // not computed. To avoid that, the Oracle can return the same value as
  // GetMinusInfinity(), thereby disabling this termination test and leaving
  // only the standard one using the Master Problem solution. However, if
  // unboundedness was to be declared when *FiBest <= *LowerBound, in this
  // case one could end up declaring the problem unbounded below. This is why
  // a value slightly smaller than *LowerBound is used instead.

  if( *FiLambda < HpINF ) {  // .. but only if Fi( Lambda ) is defined
   if( TrueLB )
    if( *FiBest - EpsLin * ABS( *FiBest ) <= *LowerBound )
     break;

   if( *FiBest <= *LowerBound * 
                  ( 1 - ( *LowerBound > 0 ? EpsLin : - EpsLin ) ) ) {
    Result = kUnbndd;
    break;
    }
   }

  // avoid the t-changing phase if Lambda1 is unfeasible- - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // note: one possible alternative t-strategy would be to set t to the
  // largest value that would have produced a feasible point, i.e.
  // t := *Alfa1 / ( - *ScPr1 )

  if( *FiLambda1 == HpINF )
   continue;
  else
   if( *FiLambda == HpINF ) {  // if reached feasibility - - - - - - - - - -
    GotoLambda1();             // go to the feasible point
    continue;                  // and start the actual minimization of Fi()
    }

  // the NS / SS decision - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  SSDone = DoSS();

  // compute the heuristic t- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  HpNum tt;
  if( ( SSDone && ( ! ( tSPar1 & 1 ) ) ) ||
      ( ( ! SSDone ) && ( tSPar1 & 2 ) ) )
   tt = Heuristic1();
  else
   tt = Heuristic2();

  if( SSDone ) {  // SS - - - - - - - - - - - - - - - - - - - - - - - - - - -
   BLOG( 1 , std::endl << " SS[" << CSSCntr << "]: DFi (" << DeltaFi
	          << ") >= m1 * Dv (" << ABS( m1 ) * Deltav << ") ~ Ht = "
	          << tt );

   tt = std::min( std::min( tMaior , t * mxIncr ) ,
		  std::max( t * mnIncr , tt ) );

   if( CSSCntr < MnSSC )  // increasing t is inhibited
    tt = t;
   else
    if( ( tSPar1 & tSP1Msk ) == kBLTTS )  // "balancing" long-term t-strategy
     if( ( DSTS <= tSPar2 * Sigma ) && ( CSSCntr < 10 ) ) {  //!! 10!
      BLOG( 1 , " ~ small D*_t( 1 )" );
      tt = t;
      }

   BLOG( 1 , std::endl );

   GotoLambda1();
   ParSS++;
   CSSCntr++;
   CNSCntr = 0;
   }
  else {        // NS - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   BLOG( 1 , std::endl << " NS[" << CNSCntr << "]: DFi (" << DeltaFi
	          << ") < m1 * Dv (" << ABS( m1 ) * Deltav
	          << ") ~ Ht = " << tt );

   tt = std::max( std::max( tMinor , t * mxDecr ) ,
		  std::min( t * mnDecr , tt ) );

   if( CNSCntr < MnNSC )  // decreasing t is inhibited
    tt = t;
   else
    if( *Alfa1 <= m3 * Sigma ) {
     BLOG( 1 , " ~ small Alfa1" );
     tt = t;
     }
    else
     switch( tSPar1 & tSP1Msk ) {
      case( kSLTTS ):
      case( kHLTTS ):
       if( vStar <= tSPar2 * EpsU * std::max( ABS( *FiLambda ) , HpNum( 1 ) )
	   ) {
	BLOG( 1 , " small v" );
	tt = t;
        }

       break;
      case( kBLTTS ):
       /*!! this version avoids problems which may occur with ill-set
	    tStar or tSPar2, but it may give worse performances with
	    "difficult" problems
       if( ( tSPar2 * DSTS >= Sigma ) && ( CNSCntr < 20 ) ) {
       !!*/

       if( tSPar2 * DSTS >= Sigma ) {
        BLOG( 1 , " ~ large D*_t( t* )" );
        tt = t;
        }
      }

   BLOG( 1 , std::endl );

   CNSCntr++;
   CSSCntr = 0;

   }   // end else( NS )- - - - - - - - - - - - - - - - - - - - - - - - - - -

  // actually update t- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
  if( tSPar1 & kEGTTS )  // endgame t-strategy: note the "/ 10"!!
   if( DSTS < EpsLin * std::max( ABS( *FiLambda ) , HpNum( 1 ) ) / 10 ) {
    tt = std::max( t * ( mxDecr + mnDecr ) / 2 , tMinor );
    BLOG( 1 , " ~ endgame, t = " << tt );
    }

  //!! the reverse should also be done: if sigma is small and D*( t* ) is
  //!! large, t should be increased --> but this would happen surely at
  //!! the beginning, it should be done only near the end

  if( ( tHasChgd = ( t != tt ) ) ) {
   CSSCntr = CNSCntr = 0;  // reset the counters as t changes
   t = tt;
   }

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  } while( ( ! MaxIter ) || ( ParIter < MaxIter ) );

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // main cycle ends here- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MaxIter && ( ParIter >= MaxIter ) && ( ! Result ) )
  Result = kStpIter;

 if( NDOt )
  NDOt->Stop();

 return( Result );

 }  // end( Bundle::Solve )

/*--------------------------------------------------------------------------*/

void Bundle::ReSetAlg( unsigned char RstLvl )
{
 if( ! ( RstLvl & RstAlg ) ) {  // reset algorithmic parameters - - - - - - -
  ParIter = ParSS = 0;     // reset iterations count
  CSSCntr = CNSCntr = 0;   // ... comprised consecutive NS/SS count
  
  if( t != tInit ) {       // reset t
   t = tInit;
   tHasChgd = true;
   }

  CmptaBPX();  // reset the dynamic number of fetched items

  // reset the precision of Fi computations
  EpsCurr = ABS( EInit );
  }

 //!! check if MPSolver != NULL !!

 if( ! ( RstLvl & RstCrr ) )  // reset the current point to all-0- - - - - - -
  SetLambda();

 if( ! ( RstLvl & ( RstSbg | RstCnt ) ) )  // reset everything- - - - - - - -
  RemoveItems();
 else
  if( ! ( RstLvl & RstSbg ) ) {  // reset subgrads (but not constrs)- - - - -
   if( Master->BSize() ) {       // if the bundle is nonempty
    if( ! Master->BCSize() )     // and it contains only subgradients
     RemoveItems();              // remove everything
    else
     for( Index i = Master->MaxName() ; i-- ; )
      if( Master->IsSubG( i ) ) {
       Oracle->Deleted( i );
       Delete( i );
       }
    }
   }
  else
   if( ! ( RstLvl & RstCnt ) ) // reset constrs (but not subgrads)- - - - - -
    if( Master->BSize() ) {    // if the bundle is nonempty
     if( Master->BSize() == Master->BCSize() )
                               // and it contains only constrs
      RemoveItems();           // remove everything
     else
      for( Index i = Master->MaxName() ; i-- ; )
       if( ! Master->IsSubG( i ) ) {
	Oracle->Deleted( i );
	Delete( i );
        }
     }

 if( ! ( RstLvl & RstFiV ) )  // reset the current value of Fi( Lambda )- - -
  *FiLambda = HpINF;

 }  // end( Bundle::ReSetAlg )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

cLMRow Bundle::ReadSol( cIndex_Set &I , Index &D )
{
 I = LamBase;
 D = LamDim;

 return( Lambda );
 }

/*--------------------------------------------------------------------------*/

cLMRow Bundle::ReadBestSol( cIndex_Set &I , Index &D )
{
 I = NULL;
 D = NumVar;

 return( KpBstL ? LmbdBst : Lambda );
 }

/*--------------------------------------------------------------------------*/

HpNum Bundle::ReadFiVal( cIndex wFi )
{
 return( WhichFi( FiLambda , wFi ) );
 }

/*--------------------------------------------------------------------------*/

HpNum Bundle::ReadBestFiVal( cIndex wFi )
{
 return( WhichFi( FiBest , wFi ) );
 }

/*--------------------------------------------------------------------------*/

cHpRow Bundle::ReadMult( cIndex_Set &I , Index &D , cIndex wFi )
{
 return( Master->ReadMult( I , D , wFi ) );
 }

/*--------------------------------------------------------------------------*/

HpNum Bundle::ReadLBMult( cIndex wFi )
{
 return( Master->ReadLBMult( wFi ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

void Bundle::AddVariables( Index NNwVrs , cLMRow IVs )
{
 if( NumVar >= MaxNumVar )          // no space for any new variable
  return;                           // return

 if( ! NNwVrs )                     // no variables to be added
  return;                           // just return

 if( NumVar + NNwVrs > MaxNumVar )  // not enough space for all the new vars
  NNwVrs = MaxNumVar - NumVar;      // put in only the first ones

 // add the variables in the Master - - - - - - - - - - - - - - - - - - - - -

 Master->AddVars( NNwVrs );

 // update Lambda[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that LamDim == NumVar if PPar2 == 0, so this assignment works for
 // both the PPar2 == 0 and the PPar2 > 0 case

 VectAssign( Lambda + LamDim , LMNum( 0 ) , NNwVrs );

 // update InctvCtr[] - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( PPar2 && PPar3 )
  VectAssign( InctvCtr + NumVar , Index( 0 ) , NNwVrs );

 if( IVs ) {  // check if any of the variables is strictly nonzero - - - - - -
  cHpNum MeD = std::max( Master->EpsilonD() , Oracle->GetBndEps() );
  Index i = 0;
  while( ( i < NNwVrs ) && ( Lambda[ i ] <= MeD ) &&
	 ( Master->IsNN( i ) || ( Lambda[ i ] >= - MeD ) ) )
   i++;

  if( i == NNwVrs )  // should have been IVs == NULL
   IVs = NULL;       // set it
  }

 if( IVs ) {  // if so, use SetLambda() to move to the new point- - - - - - -
  // note that SetLambda() reconstructs the active set (if any), making it
  // identical to the set { i : Lambda[ i ] > 0 }; thus, the new variables
  // with strictly positive initial value are put in the active set

  LMRow NewL = new LMNum[ NumVar + NNwVrs ];
  if( PPar2 )
   VectAssignB( NewL , Lambda , LamBase , NumVar , LMNum( 0 ) );
  else
   VectAssign( NewL , Lambda , NumVar );

  VectAssign( NewL + NumVar , IVs , NNwVrs );

  NumVar += NNwVrs;

  SetLambda( NewL );
  delete[] NewL;
  }
 else            // all the initial values are zero - - - - - - - - - - - - -
  if( PPar2 ) {  // add the new variables to the active set
   Index i = NumVar;
   Index_Set tLB = LamBase + LamDim;
   for( Index j = NNwVrs ; j-- ; )
    *(tLB++) = i++;

   *tLB = InINF;
   Master->AddActvSt( LamBase + LamDim , NNwVrs , LamBase );

   NumVar += NNwVrs;
   LamDim += NNwVrs;
   BHasChgd = true;
   }
  else {         // just update LamDim and NumVar
   NumVar += NNwVrs;
   LamDim = NumVar;
   }

 UpdtLowerBound();

 }  // end( Bundle::AddVariables )

/*--------------------------------------------------------------------------*/

void Bundle::RemoveVariables( cIndex_Set whch , Index hwmny )
{
 assert( ( whch && whch[ hwmny ] == InINF ) ||
	 ( ( ! whch ) && ( ! hwmny ) ) );

 // check if any of the variables is strictly nonzero - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool RmvNnz = false;
 cHpNum MeD = std::max( Master->EpsilonD() , Oracle->GetBndEps() );
 if( whch )
  if( PPar2 ) {
   cLMRow tL = Lambda;
   cIndex_Set tw = whch;
   cIndex_Set tLB = LamBase;
   for( Index h ; ( h = *(tLB++) ) < InINF ; tL++ )
    if( h == *tw ) {
     tw++;
     if( ( *tL > MeD ) || ( ( ! Master->IsNN( h ) ) && ( *tL < - MeD ) ) ) {
      RmvNnz = true;
      break;
      }
     }
   }
  else {
   cIndex_Set tw = whch;
   for( Index h ; ( h = *(tw++) ) < InINF ; )
    if( ( Lambda[ h ] > MeD ) ||
	( ( ! Master->IsNN( h ) ) && ( Lambda[ h ] < - MeD ) ) ) {
     RmvNnz = true;
     break;
     }
   }
 else {
  cLMRow tL = Lambda;
  for( Index h = 0 ; h < LamDim ; h++ , tL++ )
   if( ( *tL > MeD ) || ( ( ! Master->IsNN( h ) ) && ( *tL < - MeD ) ) ) {
    RmvNnz = true;
    break;
    }
  }

 if( RmvNnz )  // if so, first make a change of current point - - - - - - -
 {             // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
               // note that SetLambda() reconstructs the active set (if any)
               // as the set { i : Lambda[ i ] > 0 }; thus, the variables to
               // be deleted are eliminated from the active set, as required

  LMRow NewL = new LMNum[ NumVar ];

  if( whch ) {
   if( PPar2 )
    VectAssignB( NewL , Lambda , LamBase , NumVar );
   else
    VectAssign( NewL , Lambda , NumVar );

   VectAssign( NewL , LMNum( 0 ) , whch );
   }
  else
   VectAssign( NewL , LMNum( 0 ) , NumVar );

  SetLambda( NewL );
  delete[] NewL;
  }
 else          // otherwise
  if( PPar2 )    // eliminate the variables from the active set - - - - - - -
   if( LamDim ) {  // ... if the active set is nonempty - - - - - - - - - - -
    // first count how many variables will be eliminated- - - - - - - - - - -

    cIndex OldLD = LamDim;

    if( whch ) {
     cIndex_Set tw = whch;
     cIndex_Set tLB = LamBase;
     for( Index i = LamDim ; i-- ; )
      if( *(tLB++) == *tw ) {
       LamDim--;
       tw++;
       }
     }
    else
     LamDim = 0;

    // construct the new base and the set of eliminated variables - - - - - -

    if( OldLD > LamDim ) {
     Index_Set ElmV;  // the variables to be eliminated from the active set

     if( whch ) {
      ElmV = nBase + MaxNumVar + 1 - OldLD + LamDim;

      cIndex_Set tw = whch;
      Index_Set tEl = ElmV;
      Index_Set tnB = nBase;
      cIndex_Set tLB = LamBase;
      for( Index i = OldLD ; i-- ; tLB++ )
       if( *tLB == *tw )
        *(tEl++) = *(tw++);
       else
        *(tnB++) = *tLB;

      *tnB = InINF;
      }
     else {
      *nBase = InINF;
      ElmV = LamBase;
      }

     // change the active set

     Master->RmvActvSt( ElmV , OldLD - LamDim , nBase );
     Swap( LamBase , nBase );
     BHasChgd = true;
     }
    }

 // now remove the variables- - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Master->RmvVars( whch , hwmny );

 // update Lambda[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // at this point, all the variables to be eliminated have been given value
 // zero; furthermore, if the active set is used then they are all out of
 // the active set

 if( ! PPar2 ) {
  if( whch )
   Compact( Lambda , whch , NumVar );
  else
   VectAssign( Lambda , LMNum( 0 ) , NumVar );
  }

 // update NumVar - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( whch )
  NumVar -= hwmny;
 else
  NumVar = 0;

 if( ! PPar2 )
  LamDim = NumVar;

 UpdtLowerBound();

 }  // end( Bundle::RemoveVariables )

/*--------------------------------------------------------------------------*/

void Bundle::ChgFiV( cIndex wFi )
{
 HpRow tA = new HpNum[ BPar2 + 1 ];

 if( wFi > 0 ) {
  cIndex MxNm = Master->MaxName( wFi );
  if( wFi <= NrFi ) {
   if( IsEasy && IsEasy[ wFi ] )
    throw( NDOException( "Bundle::ChgFiV: wFi is easy!" ) );

   for( Index i = 0 ; i < MxNm ; i++ )
    if( Master->WComponent( i ) == wFi ) {
     cHpNum Ai = Oracle->GetVal( i );
     if( Ai == - HpINF )
      Bundle::RemoveItem( i );
     else
      tA[ i ] = Ai;
     }
   }
  else
   for( Index i = 0 ; i < MxNm ; i++ ) {
    cHpNum Ai = Oracle->GetVal( i );
    if( Ai == - HpINF )
     Bundle::RemoveItem( i );
    else
     tA[ i ] = Ai;
    }
  }

 if( ( ! wFi ) || ( wFi == InINF ) )
  tA[ BPar2 ] = Oracle->GetVal( BPar2 );

 Master->ChgAlfa( tA , wFi );

 *FiLambda = *FiBest = HpINF;

 delete[] tA;

 UpdtLowerBound();

 }  // end( Bundle::ChgFiV )

/*--------------------------------------------------------------------------*/

void Bundle::ChgSbG( cIndex strt , Index stp , cIndex wFi )
{
 Master->ChgSubG( strt , stp , wFi );

 *FiLambda = *FiBest = HpINF;

 UpdtLowerBound();

 }  // end( Bundle::ChgSbG )

/*--------------------------------------------------------------------------*/

void Bundle::RemoveItem( cIndex Name )
{
 if( OOBase[ Name ] == Inf< SIndex >() )  // no such item exists ...
  return;                                 // silently return

 Delete( Name );                          // delete it

 if( ! Master->MaxName() )                // it was the last item in the Bundle
  FreDim = 0;                             // which is now empty

 }  // end( RemoveItem )

/*--------------------------------------------------------------------------*/

void Bundle::RemoveItems( void )
{
 if( Master )
  Master->RmvItems();  // remove all items from the MPSolver (if any)

 if( Oracle )
  Oracle->Deleted();   // tell the oracle (if any) about it

 FreDim = 0;

 }  // end( RemoveItems )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

Bundle::~Bundle()
{
 // output times and statistics - - - - - - - - - - - - - - - - - - - - - - -

 BLOG( 1 , std::endl << "Total Fi() evaluations = " << FiEvaltns );

 BLOG2( 1 , NDOt , std::endl << "Tot. time (s): " << NDOt->Read() );

 if( Oracle )
  BLOG( 1 , std::endl << "Fi() time (s): " << Oracle->FiTime() );

 BLOG( 1 , std::endl );

 // memory deallocation - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Master )
  Master->SetDim();

 if( Oracle )
  MemDealloc();

 }  // end( ~Bundle )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------- HOOKS FOR DERIVED CLASSES -------------------------*/
/*--------------------------------------------------------------------------*/

Bundle::EIStatus Bundle::EveryIteration( void )
{
 return( kEINorm );

 }  // end( Bundle::EveryIteration )

/*--------------------------------------------------------------------------*/
/*----------------------- OTHER PROTECTED METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

void Bundle::FormD( void )
{
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // initialize the Master Problem Solver- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // change/set t as required- - - - - - - - - - - - - - - - - - - - - - - - -

 // Special treatment of the "empty Master Problem" case: there are no
 // subgradients, so this is only a feasibility problem, which should give
 // a feasible point as close as possible to the starting one. This is
 // "free" with some stabilizing terms (e.g. the quadratic one), but not
 // necessarily so with others (e.g. the trust region). In order to "stay
 // as close as possible", t is temporarily decreased to its minimum value.
 // As soon as there is something in the bundle, the current value of t is
 // restored (Prevt is used to hold it).

 if( Master->BCSize() >= Master->BSize() ) {
  if( ( t > tMinor ) && ( Prevt == HpINF ) ) {
   Prevt = t;
   t = tMinor;
   tHasChgd = true;
   }
  }
 else
  if( Prevt < HpINF ) {
   if( t != Prevt ) {
    t = Prevt;
    tHasChgd = true;
    }
   Prevt = HpINF;
   }

 if( tHasChgd ) {
  Master->Sett( t );
  tHasChgd = false;
  }

 // collect and set individual and global lower bounds- - - - - - - - - - - -

 if( LBHasChgd && ( *FiLambda < HpINF ) ) {
  if( *LowerBound > - HpINF ) {
   // translate the value using the reference value of the hard components,
   // but not that of the easy ones
   auto LwrBnd = *LowerBound;
      
   // the 0-th component is treated as a hard one
   auto FiHard = RfrncFi[ 0 ];
   for( Index k = 0 ; k++ < NrFi ; )
    if( ( IsEasy ) && ( IsEasy[ k ] ) )
     FiHard -= RfrncFi[ k ];
      
   LwrBnd -= FiHard;
   Master->SetLowerBound( LwrBnd );
   }
  else
   Master->SetLowerBound( - HpINF );

  for( Index k = 0 ; k++ < NrFi ; ) {
   if( IsEasy && IsEasy[ k ] )  // skip easy components
    continue;

   if( LowerBound[ k ] > - HpINF )
    Master->SetLowerBound( LowerBound[ k ] - FiLambda[ k ] , k );
   else
    Master->SetLowerBound( - HpINF , k );
  }

  LBHasChgd = false;
  }

 /* set termination criterion - - - - - - - - - - - - - - - - - - - - - - - -

 if( *FiLambda < HpINF )
  Master->SetPar( MPSolver::kZero ,
		  EpsLin * std::max( ABS( *FiLambda ) , HpNum( 1 ) )
		         / std::max( tStar / t , HpNum( 1 ) ) );
 */
          //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for(;;)  // price-in loop- - - - - - - - - - - - - - - - - - - - - - - - - -
 {        //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for(;;)  // error-handling loop - - - - - - - - - - - - - - - - - - - - - -
  {        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   // set time limit- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( MaxTime && NDOt ) {
    if( NDOt->Read() >= MaxTime ) {  // time is over already
     Result = kStpTime;
     break;
     }

    // ensure the MPSolver does not exceed the remaining time
    Master->SetPar( MPSolver::kMaxTme , MaxTime - NDOt->Read() );
    }

   MPSolver::MPStatus mps = Master->SolveMP();  // solve the MP

   if( mps == MPSolver::kOK )        // everything's alright
    break;

   if( mps == MPSolver::kUnfsbl ) {  // the feasible set is empty
    Result = kUnfsbl;
    break;
    }

   if( mps == MPSolver::kUnbndd ) {  // the MP is unbounded: this can always
                                     // be mended by decreasing t ...
    if( ( t <= tMinor ) || ( Master->BCSize() >= Master->BSize() ) ) {
     // ... but t must always be >= tMinor, and it is already == tMinor in
     // the "empty" case of the initial iteration with empty bundle
     BLOG( 1 , std::endl << "Bundle::FormD: failure in MPSolver." );
     Result = kError;
     break;
     }

    BLOG( 1 , std::endl << "Bundle::FormD: MP unbounded, decreasing t" );
    Master->Sett( t = std::max( t / 2 , tMinor ) );
    continue;
    }

   if( mps == MPSolver::kStppd ) {  // stopped by time limit
    //!! so far, the time limit in the MPSolver is only due to the
    //!! global time limit in the NDOSolver, but one day we may want
    //!! to set it independently; then, some checks will have to be
    //!! done if the solution is feasible and it can still be used
    //!! and v is sufficiently < 0: in this case we can use the
    //!! solution as well, otherwise we have to give the MPSolver
    //!! more time
    Result = kStpTime;
    break;
    }

   // mps == MPSolver::kError, i.e., there has been a numerical problem in- -
   // the MP Solver; it's not yet time to despair, as by eliminating items- -
   // it may be possible to solve the problem - - - - - - - - - - - - - - - -

   Index MBDm;
   cIndex_Set MBse;
   cHpRow Mlt = Master->ReadMult( MBse , MBDm );
   Index i = InINF;

   // the last *removable* item in Base is eliminated - - - - - - - - - - - -

   if( MBDm ) {  // if any
    if( MBse ) {
     for( ; MBDm-- ; )
      if( ( OOBase[ MBse[ MBDm ] ] >= 0 ) &&
	  ( Mlt[ MBDm ] >= Eps< HpNum >() ) ) {
       i = MBse[ MBDm ];
       break;
       }
     }
    else
     for( ; MBDm-- ; )
      if( ( OOBase[ MBDm ] >= 0 ) && ( Mlt[ MBDm ] >= Eps< HpNum >() ) ) {
       i = MBDm;
       break;
       }
    }

   if( i == InINF )  // there are no *removable* items in Base - - - -
    for( Index j = Master->MaxName() ; j-- ; )  // pick any removable item
     if( ( OOBase[ j ] >= 0 ) && ( OOBase[ j ] < Inf< SIndex >() ) ) {
      i = j;
      break;
      }

   if( i == InINF ) {  // there are no removable items at all- - - - -
    BLOG( 0 , std::endl << "Bundle::FormD: unrecoverable MP failure." );
    Result = kError;
    return;
    }

   Delete( i );      // just delete i

   }  // end ( error-handling loop )- - - - - - - - - - - - - - - - - - - - -
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Sigma = Master->ReadSigma();                  // read Sigma*

  vStar = - Master->ReadFiBLambda();            // read v*

  if( IsEasy ) {                                // there are easy component
   // adjust the contribution of the easy components to v* and Sigma*
   // easy components are treated differently from hard ones in that their
   // value in the model is *not* translated by their (reference) Fi-value in
   // Lambda; however, v* has to measure the total decrease predicted by the
   // model in Lambda1 w.r.t. the (reference) Fi-value in Lambda, and Sigma*
   // the error, again, w.r.t. the (reference) Fi-value in Lambda. since the
   // MPSolver provides non-translated Fi-values, the correction has to be
   // made here.

   HpNum EasyRifFi = 0;
   for( Index k = 0 ; k++ < NrFi ; )            
    if( IsEasy[ k ] ) {
     FiLambda1[ k ] = Master->ReadFiBLambda( k );  // read the *exact* Fi-value
     EasyRifFi += RfrncFi[ k ];
     }

   if( vStar < HpINF )
    vStar += EasyRifFi;
   Sigma += EasyRifFi;
   }

  DSTS = Master->ReadDStart( tStar );           // D_{t*,\beta,x}
  Deltav = vStar;
  if( m1 < 0 )                                  // use - z( P_{t,\beta,x} )
   Deltav -= Master->ReadDt( t );

  // Sigma* + D*_{t*}( -z* ) is the "maximum expected increase" used in
  // the stopping criterion, EpsU is that relative to Fi( Lambda )
  if( *FiLambda < HpINF )
   EpsU = ( DSTS + Sigma ) / std::max( ABS( *FiLambda ) , HpNum( 1 ) );
  else
   EpsU = 1;  // ensure EpsU is initialized somehow
  
  // the z[ i ] have changed, so in principle they are no longer in the
  // bundle: it may be the case that they actually are, but this is
  // taken care of in UpdtCntrs()
  VectAssign( whisZ + 1 , Index( InINF ) , NrFi );
  // the scalar products have changed
  VectAssign( ScPr1 , HpNum( HpINF ) , NrFi );

  if( ! PPar2 )  // no L.V.G.
   return;       // nothing else to do

  if( LamDim == NumVar )  // all variables are there
   break;                 // no "price in" to do (but possibly "price out")

  // LVG: price in- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // do pricing only for the first PPar1 iterations and then once every PPar2
  // iterations: however, do it also if convergence is detected or if the
  // subproblem is primal unfeasible

  if( ( Result == kOK ) && ( ParIter > Index( PPar1  ) ) &&
      ( ( ParIter - Index( PPar1  ) ) % PPar2 ) && ( ! IsOptimal() ) )
   return;

  // ask for *all* the entries of d[] - - - - - - - - - - - - - - - - - - - -

  cLMRow tdir = Master->Readd( true );

  // now construct the new LamBase- - - - - - - - - - - - - - - - - - - - - -

  Index nBD = 0;
  Index oBD = 0;
  LMNum epsDir = Master->EpsilonD() * t;
  Index_Set NewStuff = LamBase + LamDim;

  if( ! Master->NumNNVars() ) {  // there are no NN variables - - - - - - - -
   for( Index k = 0 ; k < NumVar ; k++ )
    if( LamBase[ oBD ] == k ) {
     oBD++;
     nBase[ nBD++ ] = k;
     }
    else
     if( ABS( tdir[ k ] ) > epsDir ) {
      nBase[ nBD++ ] = *(++NewStuff) = k;
      if( PPar3 )
       InctvCtr[ k ] = 0;

      BLOGb( LogVar , std::endl << " Created variable " << k );
      }
   }
  else
   if( Master->NumNNVars() == NumVar ) {  // there are only NN variables- - -
    for( Index k = 0 ; k < NumVar ; k++ )
     if( LamBase[ oBD ] == k ) {
      oBD++;
      nBase[ nBD++ ] = k;
      }
     else
      if( tdir[ k ] > epsDir ) {
       nBase[ nBD++ ] = *(++NewStuff) = k;
       if( PPar3 )
        InctvCtr[ k ] = 0;

       BLOGb( LogVar , std::endl << " Created variable " << k << " (>= 0)" );
       }
    }
   else {  // there are both NN and UC variables- - - - - - - - - - - - - - -
    for( Index k = 0 ; k < NumVar ; k++ )
     if( LamBase[ oBD ] == k ) {
      oBD++;
      nBase[ nBD++ ] = k;
      }
     else
      if( ( tdir[ k ] > epsDir ) ||
	  ( ( tdir[ k ] < - epsDir ) && ( ! Master->IsNN( k ) ) ) ) {
       nBase[ nBD++ ] = *(++NewStuff) = k;
       if( PPar3 )
        InctvCtr[ k ] = 0;

       BLOGb( LogVar , std::endl << " Created variable " << k );
       BLOG2b( LogVar , Master->IsNN( k ) , " (>= 0)" );
       }
    }  // end else( there are both NN and UC variables )- - - - - - - - - - -

  if( nBD == oBD )  // no changes in LamBase- - - - - - - - - - - - - - - - -
   break;
  else {  // LamBase has changed- - - - - - - - - - - - - - - - - - - - - - -
   BHasChgd = true;                // signal it
   *(++NewStuff) = InINF;          // and put termination marks to the
   nBase[ LamDim = nBD ] = InINF;  // vector of newly added stuff and
                                          // nBase
   // signal the changes to the MP solver
   Master->AddActvSt( LamBase + oBD + 1 , nBD - oBD , nBase );

   // add the entries to Lambda
   nBase--; LamBase--; Lambda--;

   for( ; nBD > oBD ; nBD-- )
    if( LamBase[ oBD ] == nBase[ nBD ] )
     Lambda[ nBD ] = Lambda[ oBD-- ];
    else
     Lambda[ nBD ] = 0;

   nBase++; LamBase++; Lambda++;

   // set the new LamBase
   Swap( LamBase , nBase );
   }

  // end price in - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  }  // end( for(;;) )

 // LVG: price out- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ParIter > Index( PPar1 ) ) {  // skip the first PPar1 iterations
  // (all of them if PPar3 == 0 ==> PPar1 == InINF )

  // ask for the entries of d[] corresponding to "active" variables - - - - -

  cLMRow tdir = Master->Readd( false );

  // now construct the new LamBase- - - - - - - - - - - - - - - - - - - - - -

  Index nBD = 0;
  Index oBD = 0;
  LMNum epsDir = Master->EpsilonD() * t;
  Index_Set DltdStuff = nBase + MaxNumVar;

  if( ! Master->NumNNVars() )  // there are no NN variables - - - - - - - - -
   for( Index k ; ( k = LamBase[ oBD ] ) < InINF ; oBD++ )
    if( ( ABS( Lambda[ oBD ] ) <= epsDir ) && ( ABS( tdir[ k ] ) <= epsDir )
	&& ( ++InctvCtr[ k ] >= Index( PPar3 ) ) ) {
     *(DltdStuff--) = k;
     BLOGb( LogVar , std::endl << " Eliminated variable " << k );
     }
    else {
     Lambda[ nBD ] = Lambda[ oBD ];
     InctvCtr[ nBase[ nBD++ ] = k ] = 0;
     }
  else
   if( Master->NumNNVars() == NumVar )  // there are only NN variables- - - -
    for( Index k ; ( k = LamBase[ oBD ] ) < InINF ; oBD++ )
     if( ( Lambda[ oBD ] <= epsDir ) && ( tdir[ k ] <= epsDir ) &&
         ( ++InctvCtr[ k ] >= Index( PPar3 ) ) ) {
      *(DltdStuff--) = k;
      BLOGb( LogVar , std::endl << " Eliminated variable " << k << " (>= 0)" );
      }
     else {
      Lambda[ nBD ] = Lambda[ oBD ];
      InctvCtr[ nBase[ nBD++ ] = k ] = 0;
      }
   else {  // there are both NN and UC variables- - - - - - - - - - - - - - -
    for( Index k ; ( k = LamBase[ oBD ] ) < InINF ; oBD++ )
     if( ( ABS( Lambda[ oBD ] ) <= epsDir )                   &&
	 ( ( tdir[ k ] <= epsDir ) &&
	   ( ( tdir[ k ] >= - epsDir ) || Master->IsNN( k ) ) ) &&
         ( ++InctvCtr[ k ] >= Index( PPar3 ) ) ) {
      *(DltdStuff--) = k;
      BLOGb( LogVar , std::endl << " Eliminated variable " << k << " (>= 0)" );
      BLOG2b( LogVar , Master->IsNN( k ) , " (>= 0)" );
      }
     else {
      Lambda[ nBD ] = Lambda[ oBD ];
      InctvCtr[ nBase[ nBD++ ] = k ] = 0;
      }
    }  // end else( there are both NN and UC variables )- - - - - - - - - - -

  if( nBD < oBD ) {  // some variables have been eliminated - - - - - - - - -
   BHasChgd = true;                // signal it
   nBase[ LamDim = nBD ] = InINF;  // and put the termination mark to nBase

   // invert the vector of deleted stuff, since it is in reverse order

   Index_Set tDSH = ++DltdStuff;
   for( Index_Set tDST = nBase + MaxNumVar ; tDSH < tDST ; )
    Swap( *(tDSH++) , *(tDST--) );

   // set the new LamBase

   Swap( LamBase , nBase );

   // signal the changes to the MP solver

   Master->RmvActvSt( DltdStuff , oBD - nBD , LamBase );
   }
  }  // end if( skip the first PPar1 iterations ) - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 }  // end( FormD )

/*--------------------------------------------------------------------------*/

void Bundle::FormLambda1( HpNum Tau )
{
 Master->MakeLambda1( Lambda , Lambda1 , Tau );

 if( ( Oracle->GetBndEps() < MPEFsb ) && Master->NumBxdVars() ) {
  // as the relative precision required to the MPSolver is not enough to
  // ensure that the bounds on the variables will be satisfied with the
  // precision required by the FiOracle, the (upper and lower) bounds are
  // strictly enforced here

  LMRow tL1 = Lambda1;

  if( Master->NumNNVars() )             // there are NN vars and UB vars
   if( Master->NumNNVars() == NumVar )  // actually, all variables are NN
    if( PPar2 ) {
     Index_Set tLB = LamBase;
     for( Index h ; ( h = *(tLB++) ) < InINF ; tL1++ ) {
      if( *tL1 < 0 )
       *tL1 = 0;

      cLMNum UBh = Oracle->GetUB( h );
      if( *tL1 > UBh )
       *tL1 = UBh;
      }
     }
    else
     for( Index h = 0 ; h < LamDim ; h++ , tL1++ ) {
      if( *tL1 < 0 )
       *tL1 = 0;

      cLMNum UBh = Oracle->GetUB( h );
      if( *tL1 > UBh )
       *tL1 = UBh;
      }
   else                                 // not all variables are NN
    if( PPar2 ) {
     Index_Set tLB = LamBase;
     for( Index h ; ( h = *(tLB++) ) < InINF ; tL1++ ) {
      if( Master->IsNN( h ) && ( *tL1 < 0 ) )
       *tL1 = 0;

      cLMNum UBh = Oracle->GetUB( h );
      if( *tL1 > UBh )
       *tL1 = UBh;
      }
     }
    else
     for( Index h = 0 ; h < LamDim ; h++ , tL1++ ) {
      if( Master->IsNN( h ) && ( *tL1 < 0 ) )
       *tL1 = 0;

      cLMNum UBh = Oracle->GetUB( h );
      if( *tL1 > UBh )
       *tL1 = UBh;
      }
  else  // there are only UB vars
   if( PPar2 ) {
    Index_Set tLB = LamBase;
    for( Index h ; ( h = *(tLB++) ) < InINF ; tL1++ ) {
     cLMNum UBh = Oracle->GetUB( h );
     if( *tL1 > UBh )
      *tL1 = UBh;
     }
    }
   else
    for( Index h = 0 ; h < LamDim ; h++ , tL1++ ) {
     cLMNum UBh = Oracle->GetUB( h );
     if( *tL1 > UBh )
      *tL1 = UBh;
     }

  }  // end( if( the bounds have to be enforced ) )

 LHasChgd = true;  // signal that Lambda1 has changed

 VectAssign( FiStatus + 1 , FiOracle::kFiStop , NrFi );

 if( BHasChgd && LamBase )
  VectAssign( Lam1Bse , LamBase , LamDim + 1 );

 VectAssign( whisG1 + 1 , Index( InINF ) , NrFi );
 // reset the representatives

 }  // end( FormLambda1 )

/*--------------------------------------------------------------------------*/

bool Bundle::FiAndGi( void )
{
 // set the new Lambda1[] and Lam1Bse[], if any - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( BHasChgd ) {  // if LamBase has changed, pass the new one to the Oracle
  if( LamDim < NumVar )
   Oracle->SetLamBase( Lam1Bse , LamDim );
  else
   Oracle->SetLamBase( NULL , NumVar );

  BHasChgd = false;
  }

 if( LHasChgd ) {  // if Lambda has changed, pass the new one to the Oracle
  Oracle->SetLambda( Lambda1 );

  LHasChgd = false;
  }

 // call Fi() - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 *FiLambda1 = Oracle->Fi( 0 );  // the 0-th component

 for( Index k = 0 ; k++ < NrFi ; )
  // now all other components, one by one - - - - - - - - - - - - - - - - - -
  if( IsEasy && IsEasy[ k ] )     // an easy component
   *FiLambda1 += FiLambda1[ k ];  // nothing to do: Fi[ k ] is known already
  else {                          // an hard component
   /* FiStatus[ k ] == kFiStop it means that the function value of component
      k still needs to be computed "with enough accuracy" (this is set in
      FormLambda1() when the function has not been computed at all); here we
      make sure to compute again only the components that have not been
      "solved accurately enough yet" (at the first call after a FormLambda1(),
      all of them). */

   if( FiStatus[ k ] == FiOracle::kFiStop ) {
    FiLambda1[ k ] = Oracle->Fi( k );
    FiStatus[ k ] = Oracle->GetFiStatus( k );
    if( FiStatus[ k ] == FiOracle::kFiError ) {
     BLOG( 1 , " ~ Error in the FiOracle" << std::endl );
     return( false );
     }
    }

   if( *FiLambda1 < HpINF ) {
    if( FiLambda1[ k ] >= HpINF )
     *FiLambda1 = HpINF;
    else
     *FiLambda1 += FiLambda1[ k ];
    }
   }

 if( *FiLambda1 == HpINF )    // Fi() is not defined in Lambda1
  DeltaFi = HpINF;
 else
  DeltaFi = *RfrncFi - *FiLambda1;

 SSDone = false;

 if( *FiLambda1 == - HpINF )  // error in computing Fi()
  return( true );

 FiEvaltns++;

 // update FiBest, if necessary - - - - - - - - - - - - - - - - - - - - - - -

 if( *FiLambda1 < *FiBest ) {
  VectAssign( FiBest , FiLambda1 , NrFi + 1 );

  if( KpBstL ) {
   if( LamDim < NumVar )
    VectAssignB( LmbdBst , Lambda1 , Lam1Bse , NumVar );
   else
    VectAssign( LmbdBst , Lambda1 , NumVar );
   }
  }

 // call (possibly many times) Gi() - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool dChanges = false;
 Index wFi = 0;  // which component of Fi the new item shall belong to; it
                 // *can* be 0, as the 0-th component (if HpINF) may
                 // generate global constraints

 for( Index Ftchd = 0 ; Ftchd < aBP3 ; ) {
  // check if the Oracle has any new items- - - - - - - - - - - - - - - - - -

  if( ! FindNextSG( wFi ) )
   break;

  // check if aggregation has to be performed - - - - - - - - - - - - - - - -
  // doing this now could occasionally result in useless aggregations, but it
  // avoids complications in the interface of MPSolver (inserting some
  // Z[ wFi ] while inserting the new item)

  Index wh = BStrategy( wFi );

  // get the space for the item from the MPSolver - - - - - - - - - - - - - -

  SgRow G1 = Master->GetItem( wFi );

  // fetch the item from the Oracle - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set SGBse;
  cIndex SGBDm = Oracle->GetGi( G1 , SGBse );

  HpNum Alfa1k = Oracle->GetVal();

  GiEvaltns++;

  // pass the base to the MP Solver - - - - - - - - - - - - - - - - - - - - -

  Master->SetItemBse( SGBse , SGBDm );

  // calculate ScPr1k and Alfa1k- - - - - - - - - - - - - - - - - - - - - - -

  Index cp;
  HpNum ScPr1k;

  if( FiLambda1[ wFi ] == HpINF )  // it is a constraint
   cp = Master->CheckCnst( Alfa1k , ScPr1k , Lambda );
  else                             // it is a subgradient
   cp = Master->CheckSubG( FiLambda1[ wFi ] - RfrncFi[ wFi ] ,
                           t , Alfa1k , ScPr1k );

  // check if the item changes the solution of the MP - - - - - - - - - - - -

  if( *FiLambda < HpINF )  // ... but only if Fi( Lambda ) is defined
   dChanges |= Master->ChangesMPSol();
  else
   dChanges = true;               // any first item changes the solution

  if( cp < InINF ) {  // the item is a copy- - - - - - - - - - - - - -
   BLOGb( LogBnd , std::endl << "New item is a copy of " << cp );

   cHpNum OrigA1k = (Master->ReadLinErr())[ cp ];

   if( OrigA1k > Alfa1k ) {        // if the copy has smaller Alfa than the
    Master->SubstItem( wh = cp );  // original, substitute it

    BLOGb( LogBnd , " with smaller Alfa" );
    }
   else
    wh = InINF;               // otherwise, nothing new has happened
   }
  else {              // insert the item, if there is space - - - - - - - - -
   if( wh < InINF )     // someone has been selected in BStrategy()
    Master->RmvItem( wh );     // remove it from the MP
   else
    wh = FindAPlace( wFi );    // find a spot in the bundle

   if( wh == InINF ) {  // no space found ...
    if( ! Ftchd ) {            // ... and this was the first item
     BLOG( 0 , std::endl << " ERROR: No space in the bundle" << std::endl );
     Result = kError;          // signal an error
     dChanges = true;          // ensure that the outer Fi-cycle ends
     }
    else
     BLOG( 1 , std::endl << " WARNING: No space in the bundle" << std::endl );

    break;                     // the cycle ends
    }

   Master->SetItem( wh );      // insert the item in the MP Solver
   OOBase[ wh ] = -1;          // ensure it won't be touched again this round

   #if( LOG_BND )
    if( NDOLLvl & LogBnd ) {
     if( FiLambda1[ wFi ] == HpINF )
      *NDOLog << std::endl << "New constraint " << wh << ", rhs = " << Alfa1k;
     else
      *NDOLog << std::endl << "New eps-subgradient " << wh << " for Fi[ "
	      << wFi << " ], eps = " << Alfa1k << ", gd = " << - ScPr1k;
     }
   #endif
   }

  // if something was inserted, bookkeeping is needed - - - - - - - - - - - -

  if( wh < InINF ) {
   Ftchd++;                    // one more item
   Oracle->SetGiName( wh );    // tell the name of the item to the FiOracle

   if( FiLambda1[ wFi ] < HpINF ) {  // it is a subgradient
    if( ( whisG1[ wFi ] == InINF ) || ( Alfa1k < Alfa1[ wFi ] ) ||
        ( ( Alfa1k == Alfa1[ wFi ] ) && ( ScPr1k > ScPr1[ wFi ] ) ) ) {
     whisG1[ wFi ] = wh;       // wh is the new representative of wFi
     Alfa1[ wFi ] = Alfa1k;
     ScPr1[ wFi ] = ScPr1k;
     }
    }
   else
    OOBase[ wh ] = -Inf< SIndex >();
    /* if the item is a constraint, mark it as permanently fixed: this may be
       a bad choice in practice, although it is required by the theory
       (we'll see ...) */
   }
  }  // end of item-collecting loop - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // compute *Alfa1 and *ScPr1 - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 *Alfa1 = 0;
 *ScPr1 = Master->ReadGid();

 for( wFi = 0 ; wFi++ < NrFi ; )
  if( whisG1[ wFi ] < InINF ) {
   if( Alfa1[ wFi ] == HpINF )
    Alfa1[ wFi ] = (Master->ReadLinErr())[ whisG1[ wFi ] ];

   *Alfa1 += Alfa1[ wFi ];

   if( ScPr1[ wFi ] == HpINF )
    ScPr1[ wFi ] = Master->ReadGid( whisG1[ wFi ] );

   *ScPr1 += ScPr1[ wFi ];
   }
  else
   Alfa1[ wFi ] = ScPr1[ wFi ] = 0;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( dChanges );

 }  // end( FiAndGi )

/*--------------------------------------------------------------------------*/

void Bundle::GotoLambda1( void )
{
 // compute the DeltaFi[] vector- - - - - - - - - - - - - - - - - - - - - - -

 VectDiff( FiLambda , FiLambda1 , RfrncFi , NrFi + 1 );

 // do the move - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Swap( Lambda , Lambda1 );
 Swap( FiLambda , FiLambda1 );
 VectAssign( RfrncFi , FiLambda , NrFi + 1 );

 // change the current point in the MP Solver - - - - - - - - - - - - - - - -

 Master->ChangeCurrPoint( t , FiLambda1 );

 // signal that Alfa1[] is not reliable - - - - - - - - - - - - - - - - - - -

 VectAssign( Alfa1 , HpNum( HpINF ) , NrFi + 1 );

 // check the new linearization errors- - - - - - - - - - - - - - - - - - - -

 CheckAlfa( true );

 // signal that the latest Fi() point is current- - - - - - - - - - - - - - -

 SSDone = true;

 }  // end( GotoLambda1 )

/*--------------------------------------------------------------------------*/

void Bundle::UpdtCntrs( void )
{
 // increase all the OOBase[] counters but those == +/-Inf< SIndex >()- - - -
 // items whose OOBase[] becomes 0 (e.g. the newly entered items, which have
 // OOBase[] == -1) are set to +1, in such a way that only the items in the
 // optimal base have OOBase[] == 0; note that the converse is not true, as
 // items in the optimal base may have OOBase[] < 0 instead

 for( SIndex_Set tOO = OOBase + Master->MaxName() ; tOO-- > OOBase ; )
  if( ( *tOO < Inf< SIndex >() ) && ( *tOO > -Inf< SIndex >() ) ) {
   (*tOO)++;
   if( ! *tOO )
    (*tOO)++;
   }

 // set to 0 the OOBase[] counter for items in base (if not < 0)- - - - - - -
 // note that there is a case in which a component wFi has Z[ wFi ] "for free"
 // in the bundle: this is when wFi only has *one* subgradient in base (or, in
 // practice, a subgradient with multiplier very close to one). This is
 // checked here (it is basically for free), and in case whisZ[] is properly
 // set so as to avoid pointless aggregations and OOBase[] is set to -1,
 // because under no circumstances such a subgradient can ever be removed
 // from the bundle

 cIndex_Set MBse;
 cHpRow Mlt = Master->ReadMult( MBse , MBDim );
 if( MBse ) {
  for( Index i ; ( i = *(MBse++) ) < InINF ; Mlt++ )
   if( *Mlt >= Eps< HpNum >() ) {
    if( ( *Mlt >= 1 - MPEOpt ) && Master->IsSubG( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ Master->WComponent( i ) ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }
  }
 else
  for( Index i = 0 ; i < MBDim ; i++ , Mlt++ )
   if( *Mlt >= Eps< HpNum >() ) {
    if( ( *Mlt >= 1 - MPEOpt ) && Master->IsSubG( i ) ) {
     // will never happen twice for the same wFi
     whisZ[ Master->WComponent( i ) ] = i;
     OOBase[ i ] = std::min( SIndex( -1 ) , OOBase[ i ] );
     }
    else
     if( OOBase[ i ] > 0 )
      OOBase[ i ] = 0;
    }

 }  // end( UpdtCntrs )

/*--------------------------------------------------------------------------*/

void Bundle::SimpleBStrat( void )
{
 for( SIndex_Set tOO = OOBase + Master->MaxName() ; tOO-- > OOBase ; )
  if( ( *tOO < Inf< SIndex >() ) && ( *tOO > SIndex( BPar1 ) ) ) {
   cIndex h = tOO - OOBase;
   Oracle->Deleted( h );
   Delete( h );
   }
 }  // end( SimpleBStrat )

/*--------------------------------------------------------------------------*/

void Bundle::Log1( void )
{
 #if( LOG_BND )
  if( NDOLLvl > 1 ) {
   *NDOLog << std::endl << "{" << SCalls << "-" << ParIter << "-"
	   << Master->MaxName() - FreDim << "-" << MBDim << "} t = " << t
	   << " ~ D*_1( z* ) = " << Master->ReadDStart( 1 )
	   << " ~ Sigma = " << Sigma << std::endl << "           ";
   if( PPar2 )
    *NDOLog << "LamDim = " << LamDim << " ~ ";

   *NDOLog <<  " Fi = ";

   if( *FiLambda == HpINF )
    *NDOLog << " - INF";
   else
    *NDOLog <<  - *FiLambda << " ~ eU = " << EpsU;

   if( BPar6 )
    *NDOLog << " ~ BP3 = " << aBP3;
   }
 #endif
 } // end( Log1 )

/*--------------------------------------------------------------------------*/

void Bundle::Log2( void )
{
 #if( LOG_BND )
  if( NDOLLvl > 1 ) {
   *NDOLog << std::endl << "           ";

   if( *LowerBound > - HpINF )
    *NDOLog << "UB = " << - *LowerBound << " ~ ";

   *NDOLog << "Fi1 = ";

   if( *FiLambda1 == - HpINF )
    *NDOLog << "+ INF => STOP." << std::endl;
   else
    if( *FiLambda1 == HpINF )
     *NDOLog << " - INF" << std::endl;
    else
     *NDOLog << - *FiLambda1 << " ~ Alfa1 = " << *Alfa1
	     << " ~ Gi1xd = " << - *ScPr1 << std::endl;
   }
 #endif
 } // end( Log2 )

/*--------------------------------------------------------------------------*/

bool Bundle::CheckAlfa( const bool All )
{
 // check for negative alfas, component-wise- - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that All == true when a Serious Step is taken, while All == false
 // when the Alfas to be checked are those just obtained in FiAndGi(); in
 // the first case, Alfa1[] is not significant and should be ignored
 // returns true is a negative Alfa was found

 Index NNegs = 0;  // number of negative Alfas found

 if( IsEasy ) {
  for( Index k = 0 ; k++ < NrFi ; )
   if( ! IsEasy[ k ] ) 
    DeltaAlfa[ k ] = - std::max( ABS( FiLambda[ k ] ) , HpNum( 1 ) ) * MPEOpt;
  }
 else
  for( Index k = 0 ; k++ < NrFi ; )
   DeltaAlfa[ k ] = - std::max( ABS( FiLambda[ k ] ) , HpNum( 1 ) ) * MPEOpt;

 if( All ) {  // check all linearization errors (after a SS)- - - - - - - - -
              //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cHpRow As = Master->ReadLinErr();
  for( Index i = 0 ; i < Master->MaxName() ; i++ , As++ ) {
   cIndex wFi = Master->WComponent( i );
   if( wFi < InINF )
    if( *As < DeltaAlfa[ wFi ] ) {
     DeltaAlfa[ wFi ] = *As;
     NNegs++;
     }
   }
  }
 else {  // only check the newly obtained stuff - - - - - - - - - - - - - - -
         // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  /* Dirty trick: we exploit the fact that in Alfa1[ wFi ] we have the
     *minimum* value of the linearization error found during the last call to
     FiAndGi(). Actually, this is not even true: Alfa1[ wFi ] is the minimum
     of *all the last calls to FiAndGi() after the last SS*. However, this is
     still the most negative Alfa of the component wFi (if any), because all
     the negative Alfas are immediately made positive, and therefore prior to
     the last call to FiAndGi() there were no (significantly) negative Alfas.
     */

  if( IsEasy ) {
   for( Index k = 0 ; k++ < NrFi ; )
    if( ( ! IsEasy[ k ] ) && ( Alfa1[ k ] < DeltaAlfa[ k ] ) ) {
     DeltaAlfa[ k ] = Alfa1[ k ];
     NNegs++;
     }
   }
  else
   for( Index k = 0 ; k++ < NrFi ; )
    if( Alfa1[ k ] < DeltaAlfa[ k ] ) {
     DeltaAlfa[ k ] = Alfa1[ k ];
     NNegs++;
     }
  }

 if( NNegs ) {  // any negative alfa has been found - - - - - - - - - - - - -
                //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  HpNum DA0 = 0;

  for( Index k = 0 ; k++ < NrFi ; ) {
   if( IsEasy && IsEasy[ k ] )
    continue;

   cHpNum DAk = DeltaAlfa[ k ];
   if( DAk < - std::max( ABS( FiLambda[ k ] ) , HpNum( 1 ) ) * MPEOpt ) {
    FiLambda[ k ] -= DAk;   // Fi( Lambda ) was incorrect
    *FiLambda -= DAk;       // ... both k and 0 (global)
    if( ! All ) {           // if it is Lambda1 to be tested
     Alfa1[ k ] -= DAk;     // correct Alfa1[] as well, both k
     *Alfa1 -= DAk;         // ... and *Alfa1
     DeltaFi -= DAk;        // DeltaFi was also wrong
     }
    DA0 -= DAk;
    DeltaAlfa[ k ] = - DAk;
    }
   else
    DeltaAlfa[ k ] = 0;
   }

  // change the Alfas in the MPSolver - - - - - - - - - - - - - - - - - - - -

  Master->ChgAlfa( DeltaAlfa );

  // FiBest and RfrncFi *may* be wrong: the only safe thing is to reset them

  VectAssign( FiBest  , FiLambda , NrFi + 1 );
  VectAssign( RfrncFi , FiLambda , NrFi + 1 );

  BLOG( 1 , " ~ Fi += " << DA0 << " (" << NNegs << ")" );

  return( true );
  }
 else
  return( false );

 }  // end( CheckAlfa )

/*--------------------------------------------------------------------------*/

void Bundle::StrongCheckAlfa( void )
{
 #if( CHECK_DS & 1 )
  if( *FiLambda < HpINF ) {  // only if Fi( Lambda ) is defined
   // set Lambda as the current point in the Oracle - - - - - - - - - - - - -

   if( ! SSDone ) {  // ... if it's not so already
    if( LamDim < NumVar )
     Oracle->SetLamBase( LamBase , LamDim );
    else
     Oracle->SetLamBase( NULL , NumVar );

    Oracle->SetLambda( Lambda );

    BHasChgd = true;
    LHasChgd = true;

    Oracle->Fi();
    }

   // now check the Alfas - - - - - - - - - - - - - - - - - - - - - - - - - -

   cHpRow MPSA = Master->ReadLinErr();

   for( Index i = 0 ; i < Master->MaxName() ; i++ ) {
    cIndex WFi = Master->WComponent( i );
    if( WFi < InINF ) {
     cHpNum Ai = Oracle->GetVal( i );
     cHpNum AEps = std::max( ABS( FiLambda[ WFi ] ) , HpNum( 1 ) ) * MPEOpt;

     if( ABS( Ai - MPSA[ i ] ) > AEps )
      if( NDOLog )
       *NDOLog << std::endl << "Error in Alfa[ " << i << " ]: MPS = " << MPSA[ i ]
	       << ", Orcl = " << Ai;
      else
       assert( false );
     }
    }
   }
 #endif

 }  // end( StrongCheckAlfa )

/*--------------------------------------------------------------------------*/

void Bundle::UpdtLowerBound( void )
{
 cHpNum LwrBnd = Oracle->GetLowerBound();
 if( LwrBnd != LowerBound[ 0 ] ) {
  LowerBound[ 0 ] = LwrBnd;
  LBHasChgd = true;
  TrueLB = ( LwrBnd > Oracle->GetMinusInfinity() );
  }

 for( Index k = 0 ; k++ < NrFi ; ) {
  if( IsEasy && IsEasy[ k ] )  // skip easy components
   continue;

  cHpNum LwrBndk = Oracle->GetLowerBound( k );
  if( LwrBndk != LowerBound[ k ] ) {
   LowerBound[ k ] = LwrBndk;
   LBHasChgd = true;
   }
  }
 } // end( UpdtLowerBound )

/*--------------------------------------------------------------------------*/

void Bundle::UpdtaBP3( void )
{
 cIndex tBP3 = ( BPar3 >= 0 ? Index( std::ceil( BPar3 ) ) :
		              Index( std::ceil( NrFi * ( - BPar3 ) ) ) );
 switch( BPar6 ) {
  case( 4 ):
   if( *FiLambda < HpINF )
    aBP3 = ( BPar5 > 0 ? aBP4 : tBP3 ) +
           Index( BPar5 / std::log10( EpsU / EpsLin ) );
    break;
  case( 3 ):
   if( *FiLambda < HpINF )
    aBP3 = ( BPar5 > 0 ? aBP4 : tBP3 ) +
                         Index( BPar5 / std::sqrt( EpsU / EpsLin ) );
   break;
  case( 2 ):
   if( *FiLambda < HpINF )
    aBP3 = ( BPar5 > 0 ? aBP4 : tBP3 ) +
           Index( BPar5 * ( EpsLin / EpsU ) );
   break;
  case( 1 ):
   if( BPar5 && ( ! ( ParIter % Index( ABS( BPar5 ) ) ) ) ) {
    if( BPar5 > 0 )
     aBP3++;
    else
     aBP3--;
    }
  }

 if( aBP3 > tBP3 )
  aBP3 = tBP3;
 else
  if( aBP3 < aBP4 )
   aBP3 = aBP4;

 }  // end( UpdtaBP3 )

/*--------------------------------------------------------------------------*/

void Bundle::CmptaBPX( void )
{
 aBP3 = ( BPar3 >= 0 ? Index( ceil( BPar3 ) ) :
                       Index( ceil( NrFi * ( - BPar3 ) ) ) );

 aBP4 = std::min( aBP3 , ( BPar4 >= 0 ? Index( ceil( BPar4 ) ) :
			                Index( ceil( NrFi * ( - BPar4 ) ) ) )
		  );

 if( BPar6 && ( BPar5 > 0 ) )
  aBP3 = aBP4;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

inline bool Bundle::DoSS( void )
{
 if( vStar == - HpINF )
  return( false );

 return( DeltaFi >= ABS( m1 ) * Deltav );
 }

/*--------------------------------------------------------------------------*/

inline void Bundle::Delete( cIndex i )
{
 if( Master ) {
  // check if this item was the "representative" for its component- - - - - -

  cIndex wFi = Master->WComponent( i );

  if( Master->IsSubG( i ) )  // it is a subgradient
   if( whisG1[ wFi ] == i )  // it is the representative of wFi
    whisG1[ wFi ] = InINF;   // a new representative is needed

  // delete the item with name `i' from the MP- - - - - - - - - - - - - - - -

  Master->RmvItem( i );
  }

 BLOGb( LogBnd , std::endl << "Item " << i << " removed" );

 // bookkeeping of internal data structures - - - - - - - - - - - - - - - - -

 HeapIns( FreList , i , FreDim++ );
 OOBase[ i ] = Inf< SIndex >();

 // compacting FreList[] if it's too big- - - - - - - - - - - - - - - - - - -
 // remove from FreList[] every name >= Master->MaxName(); note that every
 // ordered set *is* a Heap. apart from efficiency reasons, this is
 // needed because Master->MaxName() - FreDim is the only way in which the
 // Bundle can compute the number of "live" items

 cIndex MxNm = Master->MaxName();
 if( FreDim > MxNm ) {
  FreDim = 0;
  for( Index i = 0 ; i < MxNm ; i++ )
   if( OOBase[ i ] == Inf< SIndex >() )
    FreList[ FreDim++ ] = i;
  }

 assert( Master->MaxName() >= FreDim );

 }  // end( Bundle::Delete )

/*--------------------------------------------------------------------------*/

inline bool Bundle::FindNextSG( Index &wFi )
{
 // checks if the Oracle is willing to give out any subgradient of some
 // component, beginning from the current value of wFi and going round-robin

 bool GiSet = false;
 for( Index i = NrFi + 1 ; i-- ; wFi = ( wFi < NrFi ? wFi + 1 : 0 ) ) {
  if( wFi ) {
   if( IsEasy && IsEasy[ wFi ] )
    continue;
   }
  else
   if( *FiLambda1 < HpINF )    // Fi[ 0 ]( Lambda1 ) is defined
    continue;                         // Fi[ 0 ] cannot give anything new

  if( ( GiSet = Oracle->NewGi( wFi ) ) )
   break;
  }

 return( GiSet );
 }

/*--------------------------------------------------------------------------*/

inline Index Bundle::BStrategy( cIndex wFi )
{
 // this method implements the B-strategies of the code, i.e., which "old"
 // items are discarded if the bundle is full and a new item belonging to
 // component wFi has to be inserted
 // note: this is called *before* that we know if the place will actually be
 // required, so it may return InINF and nothing bad may happen

 // in particular, it returns InINF is there is plenty of space left
 // in the bundle so that no B-strategy (no removal or aggregation) is
 // required; picking a specific spot in the free space is the task of
 // FindAPlace(), which however is not called right away because the place
 // may end up not being needed
 if( FreDim || ( Master->MaxName() < Index( BPar2 ) ) )
  return( InINF );

 // there is not plenty of space- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // find the removable item with largest OOBase[]; among these with the
 // largest OOBase[], select that with largest Alfa[]
 // note: the Z[ wFi ] in the bundle (if any) are not removable and
 //       therefore cannot be selected, which in particular happens if
 //       wFi has only *one* subgradient in base

 Index wh = 0;
 cHpRow tA = Master->ReadLinErr();
 for( Index i = 0 ; ++i < Index( BPar2 ) ; )
  if( ( OOBase[ i ] > OOBase[ wh ] ) ||
      ( ( OOBase[ i ] == OOBase[ wh ] ) && ( tA[ i ] > tA[ wh ] ) ) )
   wh = i;

 if( OOBase[ wh ] < 0 )    // all items are non-removable: nothing else to
  return( InINF );  // do (except maybe complaining very loudly)
 else
  if( OOBase[ wh ] > 0 )   // a place is found
   return( wh );           // there are no problems, all done

 // wh is a basic item- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this means that *all* items of *all* components are either in base or not
 // removable, for otherwise we would have selected an item with OOBase > 0;
 // we cannot discard anything before having performed aggregation, but in
 // order to do so we also need to free some space for the Z[]

 wh = InINF;
 for( Index wFi2 = wFi ; ; ) {
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // this is done component-wise, round-robin starting from wFi in order to
  // try to keep the balance between the space allocated to each component
  if( IsEasy && IsEasy[ wFi2 ] ) {  // ... but skipping easy components
    wFi2 = ( wFi2 == NrFi ? 1 : wFi2 + 1 );
    if( wFi2 == wFi )
     break;
    else
     continue;
    }

  Index MBDm;
  cIndex_Set MBse;
  cHpRow Mlt = Master->ReadMult( MBse , MBDm , wFi2 , false );
  // note that since *all* items of *all* components are either in base or
  // not removable we can scan MBse[] for the items to be removed, possibly
  // ignoring items with Mlt[] == 0 -- but in fact not doing it because there
  // is not any

  if( whisZ[ wFi2 ] == InINF ) {
   // Z[ wFi2 ] is not already in: in principle aggregation might have to
   // be performed, so select the item to take Z[ wFi2 ] as the one with
   // smallest Mlt

   cHpRow tM = Mlt;
   Index whZ = InINF;
   HpNum tMin = HpINF;
   if( MBse ) {
    cIndex_Set tB = MBse;
    for( Index h ; ( h = *(tB++) ) < InINF ; tM++ )
     if( ( *tM < tMin ) && ( OOBase[ h ] >= 0 ) ) {
      whZ = h;
      tMin = *tM;
      }
    }
   else
    for( Index h = 0 ; h < MBDm ; h++ , tM++ )
     if( ( *tM < tMin ) && ( OOBase[ h ] >= 0 ) ) {
      whZ = h;
      tMin = *tM;
      }

   if( whZ < InINF ) {  // if there is space for Z[ wFi2 ]
    AggregateZ( Mlt , MBse , MBDm , wFi2 , whZ );  // put it in whZ
    Mlt = Master->ReadMult( MBse , MBDm , wFi2 , false );
    // and now ensure that Mlt and MBse are reliable again; this is
    // needed because AggregateZ() calls methods of the MPSolver which
    // may therefore "invalidate" temporary vectors like Mlt and MBse
    }
   else {
    // there is *no* removable item of this component: this in some sense is
    // no problem, because in this case aggregation is not needed, but on the
    // other hand it means that there is no way this component will provide
    // any place, so have to move to the next
    wFi2 = ( wFi2 == NrFi ? 1 : wFi2 + 1 );
    if( wFi2 == wFi )
     break;
    else
     continue;
    }
   }  // end if( Z[ wFi2 ] is not already in )

  // at this point, Z[ wFi2 ] is in the bundle: try to select the exiting
  // item as the one with the smallest Mult[] among these belonging to wFi2
  // this may fail, prompting to move to the next component, so one could
  // fear of having just done aggregation to no avail; yet I don't think
  // this can happen, as the only case would be that of having only one
  // removable subgradient in the wFi2 base (that is taken by whZ), but the
  // only reasonable case in which this can happen is that there is only one
  // subgradient at all, in which case it is Z[ wFi2 ] and no aggregation is
  // done (in other words, you do aggregation only if you have at least two
  // subgradients in base, and there is no reason one of them should be non
  // removable)
  HpNum tMin = HpINF;
  if( MBse ) {
   for( Index h ; ( h = *(MBse++) ) < InINF ; Mlt++ )
    if( ( *Mlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
     wh = h;
     tMin = *Mlt;
     }
   }
  else
   for( Index h = 0 ; h < MBDm ; h++ , Mlt++ )
    if( ( *Mlt < tMin ) && ( OOBase[ h ] >= 0 ) ) {
     wh = h;
     tMin = *Mlt;
     }

  if( wh < InINF )  // a place has been found
   break;                  // all done
  else {                   // all the items of wFi2 are not removable
   wFi2 = ( wFi2 == NrFi ? 1 : wFi2 + 1 );  // try the next component
   if( wFi2 == wFi )                        // if any has remained
    break;                                  // otherwise give up
   }
  }  // end for( wFi2 ) - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( wh );

 }  // end( BStrategy )

/*--------------------------------------------------------------------------*/

Index Bundle::FindAPlace( cIndex wFi )
{
 // this method is used to return the index of an available position in the
 // bundle where to store a new item belonging to "component" wFi; if there
 // are no possible positions left, then InINF is returned

 Index wh = InINF;

 if( FreDim )                               // there are deleted items
  wh = HeapDel( FreList , --FreDim );       // pick one
 else                                       // there are no deleted items ...
  if( Master->MaxName() < Index( BPar2 ) )  // ... but there is still space
   wh = Master->MaxName();                  // next name

 assert( Master->MaxName() >= FreDim );

 return( wh );

 }  // end( FindAPlace )

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::Heuristic1( void )
{
 if( *Alfa1 < Eps< HpNum >() )
  return( DeltaFi > Eps< HpNum >() ? tMaior : tMinor );
 else
  return( t * ( ( DeltaFi + *Alfa1 ) / ( 2 * *Alfa1 ) ) );
 }

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::Heuristic2( void )
{
 if( ABS( vStar - DeltaFi ) < Eps< HpNum >() )
  return( tMaior );
 else
  return( t * ( vStar / ( 2 * ( vStar - DeltaFi ) ) ) );
 }

/*--------------------------------------------------------------------------*/

inline void Bundle::InitMP( void )
{
 // this method is called only when *both* the FiOracle *and* the MPSolver
 // have been set, and it is re-called each time any one of the two changes
 // set the size- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Master->SetDim( BPar2 , Oracle , PPar2 ? true : false );

 Master->SetPar( MPSolver::kOptEps , MPEOpt );
 Master->SetPar( MPSolver::kFsbEps , MPEFsb );

 // insert the constant subgradient of the 0-th component - - - - - - - - - -

 cIndex_Set SG0Bse;
 cIndex SG0BDm = Oracle->GetGi( Master->GetItem( 0 ) , SG0Bse , BPar2 );

 Master->SetItemBse( SG0Bse , SG0BDm );
 Master->SetItem( InINF );

 tHasChgd = LBHasChgd = true;

 }  // end( InitMP )

/*--------------------------------------------------------------------------*/

inline void Bundle::AggregateZ( cHpRow Mlt , cIndex_Set MBse , Index MBDm ,
			  	cIndex wFi , cIndex whr )
{
 // note: this is *never* called in the "easy" case where aggregation is not
 // needed because there is only one subgradient in base (for the component
 // wFi) and therefore, its Mlt[] is == 1, since in this case whisZ[] has
 // been properly set in UpdtCntrs()

 // tell the Oracle what is going to happen - - - - - - - - - - - - - - - - -

 Oracle->Aggregate( Mlt , MBse , MBDm , whr );

 // ask the MPSolver the memory for keeping Z[ wFi ]- - - - - - - - - - - - -
 // note: Mlt and MBse could very well be "temporary" memory belonging to the
 // MPSolver, and any call to a method of the MPSolver may invalidate it;
 // the calls start now, and in fact MBse and Mlt are no longer used

 SgRow tZ = Master->GetItem( wFi );

 // read Z[ wFi ] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index ZBDm;
 cIndex_Set ZBse;
 Master->ReadZ( tZ , ZBse , ZBDm , wFi );

 // now pass Z[ wFi ] back to the MP Solver - - - - - - - - - - - - - - - - -

 Master->SetItemBse( ZBse , ZBDm );

 HpNum ScPri;
 HpNum Ai = Master->ReadSigma( wFi );          // its alfa is Sigma[ wFi ]

 // note that Tau == -1, meaning that Ai need not be changed since it is
 // already the linearization error in Lambda, but still the ScPri need be
 // computed
 Master->CheckSubG( 0 , -1 , Ai , ScPri );

 Master->RmvItem( whr );  // remove the old item in position whr

 Master->SetItem( whr );  // set Z[ wFi ] in position whr

 whisZ[ wFi ] = whr;      // Z[ wFi ] is in the bundle ...
 OOBase[ whr ] = -1;      // ... and it won't be removed in this iteration

 BLOGb( LogBnd , std::endl << "Aggregation performed into " << whr );

 }  // end( Bundle::AggregateZ )

/*--------------------------------------------------------------------------*/

inline HpNum Bundle::WhichFi( cHpRow FiVect , cIndex wFi )
{
 if( wFi == InINF )
  return( *FiVect );
 else
  if( ( wFi > 0 ) && ( wFi <= NrFi ) )
   return( FiVect[ wFi ] );
  else {
   cHpNum SumFi = SumV( FiVect + 1 , NrFi );
   return( wFi ? SumFi : *FiVect - SumFi );
   }
 }

/*--------------------------------------------------------------------------*/

inline void Bundle::MemDealloc( void )
{
 delete[] LmbdBst;

 delete[] ++DeltaAlfa;
 delete[] Alfa1;
 delete[] ScPr1;
 delete[] ++whisG1;

 delete[] LowerBound;
 delete[] ++FiStatus;

 delete[] FiLambda1;
 delete[] RfrncFi;
 delete[] FiBest;
 delete[] FiLambda;

 delete[] ++whisZ;
 delete[] FreList;
 delete[] OOBase;

 delete[] Lambda1;
 delete[] Lambda;

 delete[] InctvCtr;

 if( IsEasy )
  delete[] ++IsEasy;

 if( PPar2 ) {
  delete[] Lam1Bse;
  delete[] --nBase;
  delete[] --LamBase;
  }
 }  // end( MemDealloc )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Bundle.C -----------------------------*/
/*--------------------------------------------------------------------------*/
