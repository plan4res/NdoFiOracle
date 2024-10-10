/*--------------------------------------------------------------------------*/
/*--------------------------- File QPPnltMP.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/* The class QPPenaltyMP implements the Quadratic-Penalty Master Problem for
 * Bundle algorithms, using the specialized (QP) solvers implemented by the
 * [B]MinQuad class [see [B]MinQuad.h and references therein]. The class
 * conforms to the interface defined by the class MPSolver [see MPSolver.h].
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

#include "QPPnltMP.h"

#include "OPTvect.h"

#include <algorithm>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define CHECK_DS 0

/* If CHECK_DS > 0, various data structures are checked for correctness
   during the run of the algorithm, tremendously slowing down the algorithm
   but allowing to debug the thing.
   What data structures are checked is coded bit-wise in CHECK_DS:

    bit 0 (+ 1)  =>  the bounds on dir[] (if any) are checked
    bit 1 (+ 2)  =>  SubG[] and/or TSubG[] are checked
    bit 2 (+ 4)  =>  GTd[] is checked

   CHECK_DS > 0 forces asserts() to work within this unit. */

#if( CHECK_DS )
 #ifdef NDEBUG
  #undef NDEBUG
 #endif
#endif

#if( CHECK_DS )
 #include "NDOSlver.h"
#endif

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#if( ( G_IMPLM < 1 ) || ( G_IMPLM > 2 ) )

 #define HAVE_BTH 0

 // HAVE_BTH == 1 <=> both the row-wise and the column-wise representation
 // of G are available

 #define WHAT_4_Q G_IMPLM

 // if WHAT_4_Q == 0 then the column-wise representation of G is used for
 // calculating the entries of Q[][], that is the (restricted) scalar products
 // between items; if WHAT_4_Q > 0 the row-wise implementation is used instead

 #define WHAT_4_D G_IMPLM

 // if WHAT_4_D == 0 then the column-wise representation of G is used for
 // calculating the entries of z[], that is the convex combinations of entries
 // of items; if WHAT_4_D > 0 the row-wise implementation is used instead

#else

 #define HAVE_BTH 1

 #if( G_IMPLM == 1 )
  #define WHAT_4_Q 0
  #define WHAT_4_D 1
 #else
  #define WHAT_4_Q 1
  #define WHAT_4_D 0
 #endif

#endif

#if( G_IMPLM > 0 )
 #define TSGi( i ) TSubG[ i ]
#endif

#if( HV_NNVAR )
 #if( ! TWOSIDED )
  #define BxdVars NNVars
 #endif
#endif

// used in SetItem()- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#define size( x ) eDir * std::max( ABS( x ) , LMNum( 1 ) )

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static cIndex SltSize = 50;

/* In the column-wise implementation of the matrix, not all the columns are
   allocated from the beginning: rather, they are allocated in "slots" of the
   above size, i.e. SltSize subgradients are allocated/deallocated in one
   blow when needed. The same is done (even though the column-wise form of
   the matrix is not available) for the memory of MinQuad. */

#if( G_IMPLM > 0 )
 #if( G_IMPLM < 4 )
  static cIndex GTSltSize = 50;  // same as SltSize for the rows of the Bundle
 #else
  static cIndex GTSltSize = 1;   // keep memory to the bare minimum
 #endif
#endif

#if( ! HV_NNVAR )
 static cHpNum DefEpsD = 1e-6;   // value for EpsilonD()
#endif

static cIndex InINF = Inf< Index >();
static cHpNum HpINF = Inf< HpNum >();
static cLMNum LMINF = Inf< LMNum >();

/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*--  These functions are not implemented as methods of the class, since  --*/
/*--  they don't use directly its data structures.                        --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

static QuNum ScalarProductBB( cSgRow g1 , cSgRow g2 , cIndex_Set B )
{
 // equivalent to ScalarProduct( g1{B} , g2{B} , | B | )

 QuNum t = 0;
 for( Index h ; ( h = *(B++) ) < InINF ; )
  t += g1[ h ] * g2[ h ];

 return( t );
 }

/*--------------------------------------------------------------------------*/

static void VectAssignBB( LMRow g1 , cSgRow g2 , cHpNum x , cIndex_Set B )
{
 // g1{B} = x * g2{B} (element-wise)

 for( Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] = x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/

static void VectSumBB( LMRow g1 , cSgRow g2 , cHpNum x , cIndex_Set B )
{
 // g1{B} += x * g2{B} (element-wise)

 for( Index h ; ( h = *(B++) ) < InINF ; )
  g1[ h ] += x * g2[ h ];
 }

/*--------------------------------------------------------------------------*/
/*-------------------- IMPLEMENTATION OF QPPenaltyMP -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

QPPenaltyMP::QPPenaltyMP( istream *iStrm )
           : MPSolver() ,
             #if( HV_NNVAR )
	      BMinQuad() ,
             #else
	      MinQuad() ,
             #endif
	     RHSp( nullptr ) , tmpG1( nullptr ) ,
             #if( G_IMPLM > 0 )
	      TSubG( nullptr ) ,
             #endif
	     dir( nullptr )
{
 // initialize algorithmic parameters - - - - - - - - - - - - - - - - - - - -

 HpNum CtOff;
 DfltdSfInpt( iStrm , CtOff , HpNum( .01 ) );

 SetPricing( ( CtOff < 0 ? HpINF : CtOff ) );

 Index MxAdd;
 DfltdSfInpt( iStrm , MxAdd , Index( 0 ) );
 #if( HV_NNVAR )
  if( MxAdd )
   SetMaxVarAdd( MxAdd );
 #endif

 Index MxRmv;
 DfltdSfInpt( iStrm , MxRmv , Index( 0 ) );
 #if( HV_NNVAR )
  if( MxRmv )
   SetMaxVarRmv( MxRmv );
 #endif

 // set up scalar fields- - - - - - - - - - - - - - - - - - - - - - - - - - -

 tCurr = 1;
 MaxBSize = 0;
 ChkIde = false;

 #if( ! HV_NNVAR )
  EpsD = DefEpsD;
 #endif

 }  // end( QPPenaltyMP::QPPenaltyMP )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void QPPenaltyMP::SetDim( cIndex MxBSz, FiOracle *Oracle, const bool UsAvSt )
{
 // dellocate all existing data structures, if any- - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MaxBSize && ( ( MxBSz == 0 ) || Oracle ) )
  MemDealloc();

 if( MxBSz )
  if( Oracle )  // memory allocation- - - - - - - - - - - - - - - - - - - - -
  {             //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   MPSolver::SetDim( MxBSz , Oracle , UsAvSt );    // call the base class
   // this has to be done *before* the next check, as NrFi is computed there

   if( NrFi > 1 )
    throw( NDOException(
		    "QPPenaltyMP::SetDim: decomposable Fi not supported" ) );

   DimMinQuad = std::min( MaxBSize , SltSize );   // allocate [B]MinQuad
   SetMaxDim( MaxBSize , DimMinQuad , MaxSGLen );  // stuff first of all

   if( UsAvSt ) {
    ActvVrs = &InINF;
    ActvVDim = 0;
    }
   else
    ActvVrs = nullptr;

   dir = new LMNum[ MaxSGLen ];

   #if( G_IMPLM < 1 )
    TSGk = new SgNum[ MaxBSize ];
   #endif
 
   #if( G_IMPLM > 0 )
    // allocate the row-wise implementation
    TSubG = new SgRow[ MaxSGLen ];
    VectAssign( TSubG , nullptr , MaxSGLen );

    #if( HAVE_BTH )
     if( ActvVrs ) {
      Buffer = new SgRow[ MaxSGLen ];  // allocate the buffer of rows

      for( Index i = BFCntr = std::min( GTSltSize , MaxSGLen ) ; i-- ; )
       Buffer[ i ] = new SgNum[ MaxBSize ];

      NrAllRws = 0;
      }
     else
    #endif
     {
      #if( G_IMPLM > 3 )
       NrAllRws = CrrSGLen;
      #else
       // allocate the minimum necessary number of "chunks" of GTSltSize
       // rows, i.e. GTSltSize * ceil( CrrSGLen / GTSltSize )

       NrAllRws = std::min( MaxSGLen , GTSltSize * ( CrrSGLen / GTSltSize ) +
		                  ( CrrSGLen % GTSltSize ? GTSltSize : 0 ) );
      #endif

      for( Index i = 0 ; i < NrAllRws ; )
       TSubG[ i++ ] = new SgNum[ MaxBSize ];

      #if( HAVE_BTH )
       Buffer = 0;
       BFCntr = 0;
      #endif
      }
   #endif

   tmpG1 = new SgNum[ MaxSGLen ];

   #if( G_IMPLM < 3 )
   {
    SubG = new SgRow[ MaxBSize ];

    for( Index i = 0 ; i < DimMinQuad ; )
     SubG[ i++ ] = new SgNum[ MaxSGLen ];
    }
   #endif

   #if( ADD_CNST )
    ZBase = new Index[ MaxBSize + 1 ];
    ZMult = new HpNum[ MaxBSize ];
    ZBDim = 0;
   #endif

   RHSp = nullptr;
   RHSxd = HpINF;

   #if( HV_NNVAR )
    // set up variables in BMinQuad - - - - - - - - - - - - - - - - - - - - -

    NNVars = 0;
    #if( TWOSIDED )
     BxdVars = 0;
    #endif

    char csi = NNVar() | AcVar();  // by default, constrained variables are
    if( ! ActvVrs )                // at the bound; if the active set is not
     csi |= IsVar();               // used, they are also all defined

    char nsi = 0;                  // mark for unconstrained variables
    if( ! ActvVrs )
     nsi |= IsVar();

    for( Index i = 0 ; i < CrrSGLen ; i++ ) {
     const bool UCV = Oracle->GetUC( i );
     cLMNum UB = Oracle->GetUB( i );

     if( ( ! UCV ) || ( UB < LMINF ) )
      GetVars()[ i ] = csi;
     else
      GetVars()[ i ] = nsi;

     LowerBounds()[ i ] = - ( UCV ? LMINF : 0 );
     if( ! UCV )
      NNVars++;

     #if( TWOSIDED )
      if( ( ! UCV ) || ( UB < LMINF ) )
       BxdVars++;

      UpperBounds()[ i ] = UB;
     #else
      if( UB < LMINF )
       throw( NDOException(
                       "QPPenaltyMP::SetDim: upper bounds not supported" ) );
     #endif
     }

    InitialSetUp();  // have BMinQuad to "digest" the data about bounds
   #else
    // check that there are no constraints on the variables - - - - - - - - -

    for( Index i = 0 ; i < CrrSGLen ; i++ )
     if( ( ! Oracle->GetUC( i ) ) || ( Oracle->GetUB( i ) < LMINF ) )
      throw( NDOException(
	      "QPPenaltyMP::SetDim: constrained variables not supported" ) );
   #endif
   }
  else  // just change the max. bundle size and/or active set setting - - - -
  {     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   throw( NDOException(
	   "QPPenaltyMP::SetDim: changing max. bundle size not supported" ) );

   }   // end( else( Oracle ) ) - - - - - - - - - - - - - - - - - - - - - - -
       // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 else  // MxBSz == 0
  SetMaxDim( MaxBSize = 0 , 0 , 0 );  // clear all [B]MinQuad stuff

 // set up scalar fields- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CrrOrcl = Oracle;
 MinNewAlfa = HpINF;
 ZCmptd = dirCmptd = 0;

 #if( G_IMPLM >= 3 )
  Insrtd = InINF;
 #endif

 #if( ADD_CNST )
  CBZDim = InINF;
 #endif

 }  // end( QPPenaltyMP::SetDim )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::SetLowerBound( cHpNum LwBnd , cIndex wFi )
{
 if( ! wFi )
  throw( NDOException( "QPPenaltyMP::SetLowerBound: wFi == 0" ) );

 // when there is no RHS, the global and individual lower bounds have "the
 // same shape" (an all-0 subgradient) but possibly different "heigths",
 // so they are set in the same way

 if( ( wFi > 1 ) || ( ! RHSp ) ) {
  SetLB( LwBnd );
  return;
  } 

 if( LwBnd > - Inf< HpNum >() )
  throw( NDOException( "QPPenaltyMP::SetLowerBound: wFi == 1 and RHS" ) );

 }  // end( QPPenaltyMP::SetLowerBound )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

MPSolver::MPStatus QPPenaltyMP::SolveMP( void )
{
 if( MPt )
  MPt->Start();

 // reset ZCmptd, dirCmptd and RHSxz- - - - - - - - - - - - - - - - - - - - -

 RHSxd = HpINF;
 ZCmptd = dirCmptd = 0;

 // check correctness of SubG[] and TSubG[] - - - - - - - - - - - - - - - - -

 CheckSG();

 // solve the (QP)- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR && CNDVD_TMP )
  SettmpD( dir );
 #endif

 MQError qpres = SolveQP( tCurr );

 #if( HV_NNVAR && CNDVD_TMP )
  dir = GettmpD();
 #endif

 if( qpres == MinQuad::kFatal ) {
  if( MPt )
   MPt->Stop();

  return( MPSolver::kError );
  }

 if( qpres == MinQuad::kQPPrimUnbndd ) {
  if( MPt )
   MPt->Stop();

  return( MPSolver::kUnfsbl );
  }

 if( qpres == MinQuad::kStpTime ) {
  if( MPt )
   MPt->Stop();

  return( MPSolver::kStppd );
  }
 
 // reset CBZDim- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( ADD_CNST )
  CBZDim = InINF;
 #endif

 // set ZCmptd- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  #if( LAZY_D == 2 )
   if( BxdVars < CrrSGLen )
    ZCmptd = BxdVars ? 2 : 1;
   else
  #endif
    if( ActvVrs && ( ActvVDim < CrrSGLen ) )
     ZCmptd = 3;
    else
     ZCmptd = 4;
 #else
  ZCmptd = 1;
 #endif

 // special treatment for dir[] - - - - - - - - - - - - - - - - - - - - - - -
 // This part is necessary because QPPenaltyMP "cheats" on the subgradient of
 // the linear (0-th) component by directly adding it to all subgradients.
 // When the direction is computed by doing the convex combination of
 // subgradients, then the 0-th component is correctly brought in. But when
 // no subgradients are in base (that is, either the base is empty or it
 // only contains constraints) then the 0-th component has to be exogenously
 // added in; this is done now.

 if( RHSp && ( ReadBDim() <= ReadCBDim() ) ) {
  Computedir( true );
  VectSum( dir , RHSp , -tCurr , CrrSGLen );
  }

 // check correctness of GTd[]- - - - - - - - - - - - - - - - - - - - - - - -

 CheckGT();

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MPt )
  MPt->Stop();

 return( MPSolver::kOK );

 }  // end( QPPenaltyMP::SolveMP )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

HpNum QPPenaltyMP::ReadFiBLambda( cIndex wFi )
{
 HpNum vs = Readv();  // the value of the aggregated model

 if( ( wFi == InINF ) || ( vs == HpINF ) )
  return( vs );

 // compute ( - 1 / t ) * scalar product between the RHS and d*
 // however, the value of the predicted descent for the 0-th component is
 // rather "scalar product between the RHS and d*", which means that the
 // result (stored in RHSxd) has to be multiplied by ( - t )
 ComputeRHSxd();
  
 if( wFi == 0 )                   // if we only want that
  return( ( - tCurr ) * RHSxd );  // this is it

 // else: either the value of the first component alone or the aggregated
 // value excluding the 0-th component (but they are the same, there is
 // only one component), which we get by difference (- * - = +)

 return( vs + tCurr * RHSxd );
 }

/*--------------------------------------------------------------------------*/

cLMRow QPPenaltyMP::Readd( bool Fulld )
{
 // if no "active set" has been defined, Fulld == false makes no sense
 if( ( ! ActvVrs ) || ( ActvVDim == CrrSGLen ) )
  Fulld = true;

 Computedir( Fulld );

 return( dir );

 }  // end( QPPenaltyMP::Readd )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::ReadZ( LMRow tz , cIndex_Set &I , Index &D , cIndex wFi )
{
 // give it out in "dense" format - - - - - - - - - - - - - - - - - - - - - -
 
 I = 0;
 D = CrrSGLen;

 if( wFi == 0 ) {
  if( RHSp )
   VectAssign( tz , RHSp , CrrSGLen );
  else
   VectAssign( tz , LMNum( 0 ) , CrrSGLen );

  return;
  }

 // compute Z*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( ADD_CNST )
  if( CBZDim == InINF )
   if( ( CBZDim = ReadCBDim() ) )
    ComputeZBase();
  
  if( CBZDim )
   FullZ( tz , ZBase , ZMult );
  else
 #endif
  {
   CompleteZ();

   #if( HV_NNVAR )
    VectAssign( tz , SetD() , CrrSGLen );
   #else
    VectAssign( tz , dir , CrrSGLen );
   #endif
   }

 // if a RHS is set, subtract it from Z*- - - - - - - - - - - - - - - - - - -
 // ... unless one wants the "full" Z* with the RHS "inside"

 if( RHSp && ( wFi < InINF ) )
  VectSubtract( tz , RHSp , CrrSGLen );

 }  // end( QPPenaltyMP::ReadZ )

/*--------------------------------------------------------------------------*/

cHpRow QPPenaltyMP::ReadMult( cIndex_Set &I , Index &D , cIndex wFi ,
			      const bool IncldCnst )
{
 #if( ADD_CNST )
  if( ! IncldCnst ) {
   if( CBZDim == InINF )
    if( ( CBZDim = ReadCBDim() ) )
     ComputeZBase();

   if( CBZDim ) {
    I = ZBase;
    D = ZBDim;

    return( ZMult );
    }
   }
 #endif

 I = ReadBase();
 D = ReadBDim();

 return( MinQuad::ReadMult() );

 }  // end( QPPenaltyMP::ReadMult )

/*--------------------------------------------------------------------------*/

HpNum QPPenaltyMP::ReadLBMult( cIndex wFi )
{
 if( ! wFi )
  throw( NDOException( "QPPenaltyMP::ReadLBMult: wFi == 0" ) );

 // the multiplier of the global lower bound
 return( MinQuad::ReadLBMult() );

 }  // end( QPPenaltyMP::ReadLBMult )

/*--------------------------------------------------------------------------*/

HpNum QPPenaltyMP::ReadGid( cIndex Nm )
{
 CheckGT();
 ComputeRHSxd();

 return( - tCurr * ( Nm >= MaxBSize ? RHSxd : ReadGTz( Nm ) - RHSxd ) );
 }

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau )
{
 CompleteZ( false );  // ensure that everything is ready- - - - - - - - - - -

 #if( HV_NNVAR )
  if( ActvVrs && ( ActvVDim < CrrSGLen ) )
   AddDSprs( Lmbd1 , Lmbd , Tau );
  else
   AddD( Lmbd1 , Lmbd , Tau );
 #else
  if( ActvVrs && ( ActvVDim < CrrSGLen ) ) {
   cLMRow tL = Lmbd;
   LMRow tL1 = Lmbd1;
   cIndex_Set tAB = ActvVrs;
   for( Index h ; ( h = *(tAB++) ) < InINF ; )
    *(tL1++) = *(tL++) - Tau * dir[ h ];
   }
  else
   VectAdd( Lmbd1 , Lmbd , dir , LMNum( - Tau ) , CrrSGLen );
 #endif
  
 }  // end( QPPenaltyMP::MakeLambda1 )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::SensitAnals( HpNum &lp , HpNum &cp )
{
 HpNum foo;
 SensitAnals1( lp , cp , foo );
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

Index QPPenaltyMP::ReadGi( cIndex Nm , SgRow Gi , cIndex_Set &GB )
{
 if( ! IsThere( Nm ) )
  throw( NDOException( "ReadGi: item not present in the MP" ) );

 #if( G_IMPLM < 3 )
  VectAssign( Gi , SubG[ Nm ] , CrrSGLen );
 #else
  Index j = CrrSGLen;
  for( SgMat tTSG = TSubG ; j ; j-- )
   *(Gi++) = (*(tTSG++))[ Nm ];
 #endif

 GB = nullptr;
 return( CrrSGLen );

 }  // end( QPPenaltyMP::ReadGi )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

SgRow QPPenaltyMP::GetItem( cIndex wFi )
{
 #if( ADD_CNST )
  G1IsSubg = true;        // unless otherwise stated, it's a subgradient
 #endif
 G1Perz = HpINF;   // signal that G1Perz has to be computed
 #if( G_IMPLM >= 3 )
  Insrtd = InINF;  // signal that from now on the content of tmpG1[]
                          // may no longer be that of Insrtd
 #endif

 return( tmpG1 );

 }  // end( QPPenaltyMP::GetItem )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::SetItemBse( cIndex_Set SGBse , cIndex SGBDm )
{
 // "densify" the item if it is "sparse"- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( SGBse ) {
  Densify( tmpG1 , SGBse , SGBDm , CrrSGLen );

  if( ! SGBDm )  // it is an all-0 vector
   G1Perz = 0;
  }
 }  // end( QPPenaltyMP::SetItemBse )

/*--------------------------------------------------------------------------*/

Index QPPenaltyMP::CheckSubG( cHpNum DFi , cHpNum Tau , HpNum &Ai ,
			      HpNum &ScPri )
{
 // consider the RHS, if present- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( RHSp )
  VectSum( tmpG1 , RHSp , CrrSGLen );

 // check if the subgradient is identical to any other in the bundle- - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index IsIde = CheckBCopy();

 // calculate the scalar product tmpG1 * z- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ! G1Perz ) || ( Tau == 0 ) )
  // ... unless the vector is known to be all-zero, or Tau == 0 in which case
  // the scalar product is not needed
  ScPri = 0;
 else
  if( IsIde < InINF )
   ScPri = ReadGTz( IsIde );
  else {
   CompleteZ( false );

   #if( HV_NNVAR )
    ScPri = DPerG( tmpG1 );
   #else
    if( ActvVrs && ( ActvVDim < CrrSGLen ) )
     ScPri = ScalarProductBB( dir , tmpG1 , ActvVrs );
    else
     ScPri = ScalarProduct( dir , tmpG1 , CrrSGLen );
   #endif
   }

 G1Perz = ScPri;

 if( Tau != 0 ) {
  ComputeRHSxd();  // adjust the scalar product by eliminating the effect of
  ScPri -= RHSxd;  // the added RHS
  }

 // adjust Ai - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // Note that Fi( Lambda ) = RHS * Lambda + Fi[ 1 ]( Lambda ), so that
 // DeltaFi = RHS * ( Lambda1 - Lambda ) + DFi =
 //         = RHS * ( Tau / t ) * d + DFi =
 //         = - Tau * RHSxd + DFi
 // Since ScPri = G1Perz - RHSxz, the formula for Ai is
 //     Ai -= DeltaFi - G1 * ( Tau / t ) * d =
 //         = DFi - Tau * RHSxd + Tau * G1Perz = DFi + Tau * ScPri
 //
 // This is done only if Tau > 0, as Tau <= 0 signals that Ai is already the
 // linearization error in the current point

 if( Tau > 0 )
  Ai -= DFi + ScPri * Tau;

 ScPri *= - tCurr;

 if( Ai < MinNewAlfa )
  MinNewAlfa = Ai;

 Alfa1 = Ai;

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( IsIde );

 }  // end( QPPenaltyMP::CheckSubG )

/*--------------------------------------------------------------------------*/

Index QPPenaltyMP::CheckCnst( HpNum &Ai , HpNum &ScPri , cHpRow CrrPnt )
{
 #if( ADD_CNST )
  G1IsSubg = false;

  // check if the constraint is identical to any other in the bundle- - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index IsIde = CheckBCopy();

  // calculate the scalar product - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ! G1Perz )  // ... unless the vector is known to be all-zero
   ScPri = 0;
  else
   if( IsIde < InINF )
    ScPri = ReadGTz( IsIde );
   else {
    CompleteZ( false );

    #if( HV_NNVAR )
     ScPri = DPerG( tmpG1 );
    #else
     if( ActvVrs && ( ActvVDim < CrrSGLen ) )
      ScPri = ScalarProductBB( dir , tmpGi , ActvVrs );
     else
      ScPri = ScalarProduct( dir , tmpGi , CrrSGLen );
    #endif
    }

  G1Perz = ScPri;
  ScPri *= - tCurr;

  // adjust Ai- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ActvVrs && ( ActvVDim < CrrSGLen ) )
   Ai -= ScalarProduct( CrrPnt , tmpG1 , ActvVrs );
  else
   Ai -= ScalarProduct( CrrPnt , tmpG1 , CrrSGLen );

  Alfa1 = Ai;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  return( IsIde );
 #else
  throw( NDOException( "CheckCnst: adding constraints is not supported" ) );

  return( InINF );
 #endif

 }  // end( QPPenaltyMP::CheckCnst )

/*--------------------------------------------------------------------------*/

bool QPPenaltyMP::ChangesMPSol( void )
{
 #if( ADD_CNST )
  if( ! G1IsSubg ) {  // a constraint - - - - - - - - - - - - - - - - - - - -
   // this is the constraint in the MP
   //
   //  Alfa1 >= G1 * dir = - tCurr * G1 * z;
   //
   // check if it is violated by the previous optimal dir (==> dir changes)

   if( Alfa1 + tCurr * G1Perz <=
       - eR * ReadBDim() * ( ABS( tCurr * G1Perz ) + ABS( Alfa1 ) ) )
    return( false );
   }
  else              // a subgradient- - - - - - - - - - - - - - - - - - - - -
 #endif
  {
   // this subgradient corresponds to the constraint in the MP
   //
   //  v >= G1 * dir - Alfa1 = - tCurr * G1 * z - Alfa1;
   //
   // check if it is violated by the previous optimal dir (==> dir changes)
   // a subgradient with Alfa1 <= 0 always changes dir

   if( ( Alfa1 > 0 ) && ( tCurr * G1Perz + Alfa1 + Readv() >=
	  - eR * ReadBDim() * ( ABS( tCurr * G1Perz ) + Alfa1 - Readv() ) ) )
    return( false );

   }  // end( else( a subgradient ) ) - - - - - - - - - - - - - - - - - - - -

 return( true );

 }  // end( QPPenaltyMP::ChangesMPSol )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::SetItem( cIndex Nm )
{
 if( Nm == InINF ) {  // special case: the constant subgradient - - -
  ChgRHS( 0 , CrrSGLen );    // defer to the appropriate private method
  return;                    // nothing else to do
  }

 if( Nm >= DimMinQuad ) {  // if necessary, allocate new stuff- - - - - - - -
  #if( G_IMPLM < 3 )       // - - - - - - - - - - - - - - - - - - - - - - - -
   Index j = DimMinQuad;
  #endif

  do
   DimMinQuad = std::min( Index( DimMinQuad + SltSize ) , MaxBSize );
  while( Nm >= DimMinQuad );

  SetCrrBDim( DimMinQuad );

  #if( G_IMPLM < 3 )
   for( ; j < DimMinQuad ; )
    SubG[ j++ ] = new SgNum[ MaxSGLen ];
  #endif
  }

 // if it is the first subgradient ever, set Base2[]- - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 /* This used to be a good idea, but it creates a bug in the following
    situation: at the very first iteration, the all-0 direction is generated.
    Then, Fi() is called and a subgradient is found which violates some
    bound. The subgradient is inserted, and (here) a base is set out of that
    has some nonzeroes (in the components corresponding to violated bounds).
    After that, a ChangeCurrPoint() is called which is expected to do
    nothing, because d == 0, but in fact setting the base has changed d in
    the meantime, so the bounds are screwed up.

 #if( HV_NNVAR )
  #if( ADD_CNST )
   if( G1IsSubg )
  #endif
    if( ( ! MaxItemN() ) && NNVars ) {  // the first subgradient ever - - - -
     SgRow tGi = tmpG1;
     char *tV = GetVars();
     cLMRow tLB = LowerBounds();
     cHpNum eDir = EpsilonD() * max( BDim , Index( 1 ) );
     #if( TWOSIDED )
      cLMRow tUB = UpperBounds();
      for( Index cnt = NNVars ; cnt ; tV++ , tLB++ , tUB++ , tGi++ )
     #else
      for( Index cnt = NNVars ; cnt ; tV++ , tLB++ , tGi++ )
     #endif
       if( *tV & NNVar() ) {
        cnt--;
        if( *tV & IsVar() ) {
         cLMNum Gii = - tCurr * (*tGi);
         if( Gii + size( Gii ) < *tLB )
          *tV |= AcVar();
         else
          #if( TWOSIDED )
	   if( Gii - size( Gii ) > *tUB )
	    *tV |= ( AcVar() | UBVar() );
	   else
          #endif
	    *tV &= ~AcVar();
        }
       }

     InitialSetUp();
     }
 #endif
 */

 // update the row-wise data structure, if necessary- - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( G_IMPLM > 0 )
  SgRow tGi = tmpG1;
  #if( HAVE_BTH )
   if( ActvVrs && ( ActvVDim < CrrSGLen ) ) {
    cIndex_Set AVt = ActvVrs;
    for( Index h ; ( h = *(AVt++) ) < InINF ; )
     TSubG[ h ][ Nm ] = tGi[ h ];
    }
   else
  #endif
   {
    SgMat tTSG = TSubG + CrrSGLen;
    for( tGi += CrrSGLen ; tTSG > TSubG ; )
     (*(--tTSG))[ Nm ] = *(--tGi);
    }
 #endif

 // update the column-wise data structure, if necessary - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( G_IMPLM < 3 )
  Swap( tmpG1 , SubG[ Nm ] );
 #endif

 // pass the item to the (QP) solver- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( ADD_CNST )
  if( ! G1IsSubg )
   AddConstr( Nm , Alfa1 );
  else
 #endif
   AddSubGrad( Nm , Alfa1 );

 SetGTz( Nm , G1Perz );  // pass the scalar product to the QP solver
 #if( G_IMPLM >= 3 )
  Insrtd = Nm;
 #endif

 }  // end( QPPenaltyMP::SetItem )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::SubstItem( cIndex Nm )
{
 if( Nm >= DimMinQuad )
  throw( NDOException( "QPPenaltyMP::SubstItem: invalid Nm" ) );

 ChangeAlfa( Nm , Alfa1 - ReadAlfa()[ Nm ] );
 #if( G_IMPLM >= 3 )
  Insrtd = Nm;
 #endif

 }  // end( QPPenaltyMP::SubstItem )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::RmvItem( cIndex i ) { ResetBundle( i ); }

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::RmvItems( void ) { ResetBundle(); }

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::SetActvSt( cIndex_Set AVrs , cIndex AVDm )
{
 if( ! ActvVrs )
  throw( NDOException( "SetActvSt: active set not initialized" ) );

 // first, take care of ActvVrs - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( AVrs ) {
  ActvVrs = AVrs;
  ActvVDim = AVDm;
  }
 else {
  ActvVrs = &InINF;
  ActvVDim = 0;
  }

 // now set up the data structures in the (QP) solver - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  char *vs = GetVars();
  cIndex_Set Actv = ActvVrs;
  for( Index i = 0 ; i < CrrSGLen ; i++ )
   if( *Actv == i ) {
    Actv++;
    *(vs++) |= IsVar();  // the variable is in (whichever type it is)

    #if( HAVE_BTH )
     if( ! TSubG[ i ] ) {
      GetNewGTRow( i );
      SetGTRow( i );
      }
    #endif
    }
   else {
    *(vs++) &= NNVar();  // the variable is out (but keep its type)

    #if( HAVE_BTH )
     if( TSubG[ i ] )
      PutGTRow( i );
    #endif
    }

  InitialSetUp();
 #else
  #if( HAVE_BTH )
   cIndex_Set Actv = ActvVrs;
   for( Index i = 0 ; i < CrrSGLen ; i++ )
    if( *Actv == i ) {
     Actv++;

     if( ! TSubG[ i ] ) {
      GetNewGTRow( i );
      SetGTRow( i );
      }
     }
    else
     if( TSubG[ i ] )
      PutGTRow( i );
  #endif

  ChangeQ();  // update Q[][] ...
  UpdateB();  // ... and the base in response to the change
 #endif

 RHSxd = HpINF;

 }  // end( QPPenaltyMP::SetActvSt )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::AddActvSt( cIndex_Set Addd , cIndex AdDm , cIndex_Set AVrs )
{
 if( ! ActvVrs )
  throw( NDOException( "AddActvSt: active set not initialized" ) );

 // add the variables to the active set - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ActvVrs = AVrs;
 ActvVDim += AdDm;

 // now add all the variables to the subproblem - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if necessary, allocate the memory for the new rows- - - - - - - - - - - -

 #if( HAVE_BTH )
  cIndex_Set tA = Addd;
  for( Index i ; ( i = *(tA++) ) < InINF ; ) {
   GetNewGTRow( i );
   SetGTRow( i );
   }
 #endif

 // add the variables to the subproblem - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  BMinQuad::AddVars( Addd , AdDm );
 #else
  const bool UptdB = ( AdDm < 2 * ReadBDim() / 3 );  // if updating the base
                                                     // is convenient
  for( Index h ; ( h = *(Addd++) ) < InINF ; )
   AddSGSpaceDim( TSGi( h ) , UptdB );

  if( ! UptdB )  // the base has not been updated
   UpdateB();    // do it now
 #endif

 RHSxd = HpINF;

 }  // end( QPPenaltyMP::AddActvSt )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::RmvActvSt( cIndex_Set Rmvd , cIndex RmDm , cIndex_Set AVrs )
{
 if( ! ActvVrs )
  throw( NDOException( "RmvActvSt: active set not initialized" ) );

 // remove the variables from the active set- - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ActvVrs = AVrs;
 ActvVDim -= RmDm;

 // now remove all the variables from the subproblem- - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // remove the variables from the subproblem- - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  BMinQuad::RemoveVars( Rmvd , RmDm );
 #else
  if( RmDm < ActvVDim ) {  // updating Q is convenient
   const bool UptdB = ( RmDm < 2 * ReadBDim() / 3 );  // if updating the base
                                                      // is also convenient
   cIndex_Set tR = Rmvd;
   for( Index h ; ( h = *(tR++) ) < InINF ; )
    CutSGSpaceDim( TSGi( h ) , UptdB );

   if( ! UptdB )  // the base has not been updated
    UpdateB();    // do it now
   }
  else {                   // recomputing Q from scratch is convenient
   ChangeQ();                       // ... so recompute it
   if( ActvVDim + 2 < ReadBDim() )  // a clear sign that the old base is bad
    ResetB();                       // restart from scratch
   else
    UpdateB();                      // update the base
   }
 #endif

 // give back the memory- - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HAVE_BTH )
  for( ; *Rmvd < InINF ; Rmvd++ )
   PutGTRow( *Rmvd );
 #endif

 RHSxd = HpINF;

 }  // end( QPPenaltyMP::RmvActvSt )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::AddVars( cIndex NNwVrs )
{
 if( ! NNwVrs )  // adding 0 new variables
  return;        // all is done

 cIndex NewSGLen = CrrSGLen + NNwVrs;
 if( NewSGLen > MaxSGLen )
  throw( NDOException( "QPPenaltyMP::AddVars: too many variables" ) );

 // if necessary, allocate the memory for the new rows- - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that, if both the row-wise & column-wise implementation of G are
 // available (HAVE_BTH > 0), then the new rows are allocated only if the
 // active set is *not* defined; otherwise, memory for the new rows will
 // be allocated only if and when they are inserted in the active set

 #if( G_IMPLM > 0 )
  #if( HAVE_BTH )
   if( ! ActvVrs )
  #endif
    while( NewSGLen > NrAllRws )
     #if( G_IMPLM <= 3 )
      for( Index j = std::min( MaxSGLen - NrAllRws , GTSltSize ) ; j-- ; )
     #endif
       TSubG[ NrAllRws++ ] = new SgNum[ MaxBSize ];
 #endif

 // ask the FiOracle the new entries of the RHS - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( G_IMPLM >= 3 )
  Insrtd = InINF;  // since tmpG1[] is used as a temporary, signal
                   // that it cannot be used any longer as G[ Insrtd ]
 #endif

 if( RHSp ) {  // if a RHS is set already - - - - - - - - - - - - - - - - - -
  cIndex_Set SGBse;
  Index SGBDim = CrrOrcl->GetGi( RHSp + CrrSGLen , SGBse , MaxBSize ,
				 CrrSGLen , NewSGLen );
  if( SGBse )
   Densify( RHSp , SGBse , SGBDim , NewSGLen , CrrSGLen );

  RHSxd = HpINF;
  }
 else {      // a non-0 RHS may have to be set- - - - - - - - - - - - - - - -
  cIndex_Set SGBse;
  Index SGBDim = CrrOrcl->GetGi( tmpG1 , SGBse , MaxBSize ,
				 CrrSGLen , NewSGLen );

  if( SGBDim ) {  // the RHS becomes nonzero now
   RHSp = new SgNum[ MaxSGLen ];
   VectAssign( RHSp , SgNum( 0 ) , CrrSGLen );
   if( SGBse ) {
    cSgRow ttG1 = tmpG1;
    for( Index i = CrrSGLen ; i < NewSGLen ; i++ )
     if( *SGBse == i ) {
      RHSp[ i ] = *(ttG1++);
      SGBse++;
      }
     else
      RHSp[ i ] = 0;
    }
   else
    VectAssign( RHSp + CrrSGLen , tmpG1 , NewSGLen - CrrSGLen );

   RHSxd = HpINF;
   }
  }

 // for each item in the bundle, ask the FiOracle the new information - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < MaxItemN() ; i++ ) {
  if( ! IsThere( i ) )  // only ask *existing" items
   continue;

  cIndex_Set SGBse;
  #if( G_IMPLM < 3 )
   SgRow tGi = SubG[ i ];
  #else
   SgRow tGi = tmpG1;
  #endif

  // ask the components [CrrSGLen , NewSGLen) for item i- - - - - - - - - - -

  Index SGBDim = CrrOrcl->GetGi( tGi + CrrSGLen , SGBse , i , CrrSGLen ,
				 NewSGLen );

  // if they are given in sparse format, densify them - - - - - - - - - - - -

  if( SGBse )
   Densify( tGi , SGBse , SGBDim , NewSGLen , CrrSGLen );

  tGi += CrrSGLen;

  // add the RHS (if any) - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( RHSp )
   VectSum( tGi , RHSp + CrrSGLen , NewSGLen - CrrSGLen );

  // copy the information in the row-wise description - - - - - - - - - - - -

  #if( G_IMPLM > 0 )
   #if( HAVE_BTH )
    if( ! ActvVrs )  // the new variables are all *not* active
   #endif
     for( Index j = CrrSGLen ; j < NewSGLen ; )
      TSubG[ j++ ][ i ] = *(tGi++);
  #endif

  }  // end( for( i ) )

 // set the sign & bounds of the variables- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index j = CrrSGLen ; j < NewSGLen ; j++ ) {
  const bool UCV = CrrOrcl->GetUC( j );
  cLMNum UB = CrrOrcl->GetUB( j );

  #if( HV_NNVAR )  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   GetVars()[ j ] = ( ( ! UCV ) || ( UB < LMINF ) ? NNVar() : 0 );
   LowerBounds()[ j ] = - ( UCV ? LMINF : 0 );
   if( ! UCV )
    NNVars++;

   #if( TWOSIDED )
    if( ( ! UCV ) || ( UB < LMINF ) )
     BxdVars++;

    UpperBounds()[ j ] = UB;
   #else
    if( UB < LMINF )
     throw( NDOException( "AddVars: upper bounds are not supported" ) );
   #endif
  #else            // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   if( ( ! UCV ) || ( UB < LMINF ) )
    throw( NDOException( "AddVars: boxed variables are not supported" ) );

   if( ! ActvVrs )               // if the active set is not used
    AddSGSpaceDim( TSGi( j ) );  // add the variable to the MP
  #endif

  }  // end( for( CrrSGLen ) )

 #if( HV_NNVAR )
  if( ! ActvVrs )                           // if the active set is not used
   BMinQuad::AddVars( CrrSGLen , NNwVrs );  // add the variables to the MP
 #endif

 // update CrrSGLen - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CrrSGLen = NewSGLen;

 }  // end( QPPenaltyMP::AddVars )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::RmvVars( cIndex_Set whch , Index hwmny )
{
 if( ! whch )  // remove all the variables- - - - - - - - - - - - - - - - - -
 {             // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ! ActvVrs ) {  // if the active set is not used, delete the variables
   #if( HV_NNVAR )
    for( char *vs = GetVars() ; CrrSGLen ; CrrSGLen-- )
     *(vs++) &= NNVar();

    InitialSetUp();
    BxdVars = NNVars = 0;
   #else
    CrrSGLen = 0;
    ChangeQ();
    ResetB();
   #endif
   }

  return;
  }

 if( ! hwmny )  // actually, no variables have to be removed
  return;       // all is done

 // update NNVars and BxdVars - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
 {
  cIndex_Set tw = whch;
  for( Index h ; ( h = *(tw++) ) < InINF ; )
   if( GetVars()[ h ] & NNVar() ) {
    #if( TWOSIDED )
     BxdVars--;

     if( LowerBounds()[ h ] > - LMINF )
    #endif
      NNVars--;
    }
  }
 #endif

 // if the active set is not used, delete the variables - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // otherwise, RmvActvSt() has to be called before for those variables that
 // happen to be in the active set (but no check is done)

 if( ! ActvVrs ) {
  #if( HV_NNVAR )
   BMinQuad::RemoveVars( whch , hwmny );
  #else
   cIndex newCSGL = CrrSGLen - hwmny;
   if( hwmny < newCSGL ) {  // updating Q is convenient
    const bool UptdB = ( hwmny < 2 * ReadBDim() / 3 );  // if updating the
                                                        // base is convenient
    cIndex_Set tw = whch;
    for( Index h ; ( h = *(tw++) ) < InINF ; )
     CutSGSpaceDim( TSGi( h ) , UptdB );

    if( ! UptdB )  // the base has not been updated
     UpdateB();    // do it now
    }
   else {                   // recomputing Q from scratch is convenient
    ChangeQ();                      // ... so recompute it

    if( newCSGL + 2 < ReadBDim() )  // a clear sign that the old base is bad
     ResetB();                      // restart from scratch
    else
     UpdateB();                     // update the base
    }
  #endif
  }

 // update the RHS- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( RHSp )
  Compact( RHSp , whch , CrrSGLen );

 // update the bundle (rows)- - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( G_IMPLM > 0 )
 {
  cIndex_Set tw = whch;

  #if( HAVE_BTH )
   if( ActvVrs ) {
    Index i = *tw;
    for( Index j = i ; j < CrrSGLen ; j++ )
     if( *tw == j ) {
      tw++;
      if( TSubG[ j ] ) {  // give back the row to the Buffer
       Buffer[ BFCntr++ ] = TSubG[ j ];
       TSubG[ j ] = 0;
       NrAllRws--;
       }
      }
     else
      TSubG[ i++ ] = TSubG[ j ];
    }
   else
  #endif
   {
    #if( G_IMPLM >= 4 )
     for( Index i ; ( i = *(tw++) ) < InINF ; )
      delete[] TSubG[ i ];

     NrAllRws -= hwmny;
     Compact( TSubG , whch , CrrSGLen );
     VectAssign( TSubG + NrAllRws , SgRow( 0 ) , hwmny );
    #else
     Index i = *(tw++);
     for( Index j = i ; ++j < CrrSGLen ; )
      if( *tw == j )
       tw++;
      else
       Swap( TSubG[ i++ ] , TSubG[ j ] );
    #endif
    }
  }
 #endif

 // update the bundle (columns) - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( G_IMPLM < 3 )
  for( Index h = 0 ; h < MaxItemN() ; h++ )
   if( IsThere( h ) )
    Compact( SubG[ h ] , whch , CrrSGLen );
 #else
  if( Insrtd < InINF )
   Compact( tmpG1 , whch , CrrSGLen );
 #endif

 // update CrrSGLen - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 CrrSGLen -= hwmny;

 RHSxd = HpINF;
 
 // rename the variables in [B]MinQuad- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // this has to be done *after* updating the Bundle for the case where the
 // data structures are checked within RenameVars()

 #if( HV_NNVAR )
  RenameVars( whch );
 #else
  Compact( dir , whch , CrrSGLen );
 #endif

 }  // end( QPPenaltyMP::RmvVars )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::ChgAlfa( cHpRow DeltaAlfa )
{
 ChangeAlfa( DeltaAlfa[ 1 ] );

 }  // end( QPPenaltyMP::ChgAlfa( cHpRow ) )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::ChgAlfa( cHpRow NewAlfa , cIndex wFi )
{
 cHpRow tA = ReadAlfa();
 for( Index i = 0 ; i < MaxItemN() ; i++ )
  if( IsThere( i ) )
   ChangeAlfa( i , NewAlfa[ i ] - tA[ i ] );

 }  // end( QPPenaltyMP::ChgAlfa( cHpRow , cIndex ) )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::ChgAlfa( cIndex i , cHpNum Ai )
{
 if( ( i >= MaxItemN() ) || ( ! IsThere( i ) ) )
  throw( NDOException( "invalid index in ChgAlfa( i )" ) );

 ChangeAlfa( i , Ai - ReadAlfa()[ i ] );

 }  // end( QPPenaltyMP::ChgAlfa( cIndex , cHpNum ) )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::ChangeCurrPoint( cLMRow DLambda , cHpRow DFi )
{
 // first, update the bounds - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  LMRow lb = LowerBounds() + CrrSGLen;
  for( cLMRow tDL = DLambda + CrrSGLen ; tDL-- > DLambda ; )
   if( *(--lb) > - LMINF )
    *lb -= *tDL;

  #if( TWOSIDED )
   LMRow ub = UpperBounds() + CrrSGLen;
   for( cLMRow tDL = DLambda + CrrSGLen ; tDL-- > DLambda ; )
    if( *(--ub) < LMINF )
     *ub -= *tDL;
  #endif

  ChangeBounds();
 #endif

 // check if the bounds are correct- - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( HV_NNVAR )
  CheckBounds();
 #endif

 // now change the Alfa of each item in the Bundle- - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // Important note: QPPenaltyMP "cheats" somewhat w.r.t. the MPSolver
 // interface, in that the latter requires to keep, together with a
 // subgradient g[ i ], its linearization error alfa[ i ] *with respect to
 // the component of Fi() that g[ i ] belongs to*. That is, if g[ i ] belongs
 // to the h-th component, then
 //
 //  alfa[ i ] = Fi[ h ]( lambda ) - Fi[ h ]( lambda_i ) -
 //              g[ i ]( lambda - lambda_i )
 //
 // where lambda is the "current" point and lambda_i is the point where
 // g[ i ] had been obtained. However, QPPenaltyMP actually keeps the
 // linearization errors *of ( g[ i ] + RHS )* of the *full function*
 // Fi[ 0 ] + Fi[ 1 ], where Fi[ 0 ]( lambda ) = RHS * lambda is the linear
 // 0-th component, that is
 //
 //  my_alfa[ i ] = ( Fi[ 0 ] + Fi[ 1 ] )( lambda ) -
 //                 ( Fi[ 0 ] + Fi[ 1 ] )( lambda_i ) -
 //                 ( RHS + g[ i ] )( lambda - lambda_i )
 //  = alfa[ i ] + Fi[ 0 ]( lambda ) - Fi[ 0 ]( lambda_i )
 //              - RHS( lambda - lambda_i )
 //
 // (note that there can be only one nonlinear component). Of course, the
 // two numbers are usually identical - the only (sub)gradient of a linear
 // function can only have zero linearization error. However, an exception
 // exists: when the actual value of Fi[ 0 ] + Fi[ 1 ] is not known, an
 // estimate is used for Fi[ 0 ]( lambda ), which may happen to be wrong.
 //
 // Thus, when ChangeCurrPoint() is called it is necessary to check that
 //
 //  alfa_0 = Fi[ 0 ]( < new point > ) - Fi[ 0 ]( < old point > )
 //           - RHS * DLambda == 0
 //
 // and, if not, update all alfa[ i ] (i.e., my_alfa[ i ]) by subtracting
 // them alfa_0, i.e., turn them back to the "true" alfa[ i ].
 //
 // Note that the value of Fi[ 0 ] is not directly available in DFi[]; in
 // fact, DFi[ 0 ] =  ( Fi[ 0 ] + Fi[ 1 ] )( < new point > ) -
 //                   ( Fi[ 0 ] + Fi[ 1 ] )( < old point > )
 // while DFi[ 1 ] =  Fi[ 1 ]( < new point > ) - Fi[ 1 ]( < old point > )
 //
 // As a consequence,
 //
 // Fi[ 0 ]( < new point > ) - Fi[ 0 ]( < old point > ) = DFi[ 0 ] - DFi[ 1 ]
  
 if( MaxItemN() ) {  // ... if any
  HpNum tDFi = *DFi;

  if( RHSp ) {
   HpNum alfa0;
   if( ActvVrs && ( ActvVDim < CrrSGLen ) )
    alfa0 = ScalarProductBB( DLambda , RHSp , ActvVrs );
   else
    alfa0 = ScalarProduct( DLambda , RHSp , CrrSGLen );

   alfa0 -= DFi[ 0 ] - DFi[ 1 ];

   if( ABS( alfa0 ) > Eps< HpNum >() * CrrSGLen * 100 ) {
    tDFi += alfa0;

    if( MPLLvl > 1 )
     *MPLog << "Warning: alfa0 = " << alfa0 << " [DLambda]" << endl;
    }
   }

  #if( G_IMPLM >= 3 )   /*- - - - - - - - - - - - - - - - - - - - - - - - -*/

  for( Index i = 0 ; i < MaxItemN() ; i++ )
   if( IsASubG( i ) )
    tmpa[ i ] = tDFi;
   else
    tmpa[ i ] = 0;

   for( Index i = 0 ; i < CrrSGLen ; i++ )
    if( DLambda[ i ] )
     VectSum( tmpa , TSubG[ i ] , - DLambda[ i ] , MaxItemN() );

   for( Index i = 0 ; i < MaxItemN() ; i++ )
    if( IsThere( i ) )
     ChangeAlfa( i , tmpa[ i ] );

  #else   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   for( Index i = 0 ; i < MaxItemN() ; i++ )
    if( IsThere( i ) )
     ChangeAlfa( i , ( IsASubG( i ) ? tDFi : 0 ) -
		     ScalarProduct( DLambda , SubG[ i ] , CrrSGLen ) );

  #endif   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  }  // end if( MaxItemN() )

 MinNewAlfa = HpINF;
 // reset the min negative Alfa, since they have changed

 }  // end( ChangeCurrPoint( DLambda , DFi ) )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::ChangeCurrPoint( cHpNum Tau , cHpRow DFi )
{
 #if( RECOMPUTE_GTZ > 0 )
  // to avoid accumulation of numerical errors in the linearization errors,
  // recompute the scalar products G_i * z* from scratch now, except if
  // Tau == 0 since then the scalar products are not used anyway

  if( Tau > 0 ) {
   #if( ! HV_NNVAR )
    cLMRow zz = dir;
   #else
    cLMRow zz = SetD();
   #endif

   #if( G_IMPLM >= 3 )  /*- - - - - - - - - - - - - - - - - - - - - - - - -*/

    VectAssign( tmpa , 0 , MaxItemN() );

    for( Index i = 0 ; i < CrrSGLen ; i++ )
     if( zz[ i ] )
      VectSum( tmpa , TSubG[ i ] , zz[ i ] , MaxItemN() );

    for( Index i = 0 ; i < MaxItemN() ; i++ )
     if( IsThere( i ) )
      SetGTz( i , tmpa[ i ] );

   #else  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

    for( Index i = 0 ; i < MaxItemN() ; i++ )
     if( IsThere( i ) )
      SetGTz( i , ScalarProduct( zz , SubG[ i ] , CrrSGLen ) );

   #endif  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   }
 #endif
 
 CheckGT();

 // change the linearization errors

 MoveAlongD( Tau , *DFi );

 // check that the linearization error of the 0-th component is zero, and
 // if it is not correct all the linearization errors; see the comments
 // within ChangeCurrPoint( DLambda , DFi ) above

 ComputeRHSxd();
 HpNum alfa0 = - tCurr * RHSxd - ( DFi[ 0 ] - DFi[ 1 ] );

 if( ABS( alfa0 ) > Eps< HpNum >() * CrrSGLen * 100 ) {
  ChangeAlfa( alfa0 );

  if( MPLLvl > 1 )
   *MPLog << "Warning: alfa0 = " << alfa0 << " [Tau]" << endl;
  }

 #if( HV_NNVAR )
  CheckBounds();
 #endif

 MinNewAlfa = HpINF; 
 // reset the min negative Alfa, since they have changed

 }  // end( QPPenaltyMP::ChangeCurrPoint( Tau , DFi ) )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::ChgSubG( cIndex strt , Index stp , cIndex wFi )
{
 if( stp > CrrSGLen )
  stp = CrrSGLen;

 if( strt >= stp )
  return;

 if( ! wFi ) {  // special case: changing the RHS only- - - - - - - - - - - -
  // recover the new entries of the RHS
  cIndex_Set SGBse;
  Index SGBDim = CrrOrcl->GetGi( tmpG1 , SGBse , MaxBSize , strt , stp );

  if( ! SGBDim )           // an all-zero vector
   G1Perz = 0;             // let this be known
  else {                   // some entries are nonzero
   G1Perz = HpINF;  // let this be known
   if( SGBse )             // if the vector is still sparse, densify it
    Densify( tmpG1 - strt , SGBse , SGBDim , stp , strt );
    // note  ^^^^^^^^^^^^ this: because there is   ^^^^ this Densify will
    // not touch the first strt entries of the vector it is given, hence
    // it will densify in place only the first stp - strt entries of
    // tmpG1 despite the indices in SGBse being (necessarily) >= strt
   }

  ChgRHS( strt , stp );    // defer to the appropriate private method
  return;                  // nothing else to do
  }

 // ask the FiOracle the new entries of the RHS - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( G_IMPLM >= 3 )
  Insrtd = InINF;  // since tmpG1[] is used as a temporary, signal
                   // that it cannot be used any longer as G[ Insrtd ]
 #endif

 if( wFi == InINF ) {
  if( RHSp ) {  // if a RHS is set already- - - - - - - - - - - - - - - - - -
   cIndex_Set SGBse;
   Index SGBDim = CrrOrcl->GetGi( RHSp + strt , SGBse , MaxBSize ,
				  strt , stp );

   if( SGBse )
    Densify( RHSp , SGBse , SGBDim , stp , strt );

   RHSxd = HpINF;
   }
  else {      // a non-0 RHS may have to be set - - - - - - - - - - - - - - -
   cIndex_Set SGBse;
   Index SGBDim = CrrOrcl->GetGi( tmpG1 , SGBse , MaxBSize , strt , stp );

   if( SGBDim ) {  // the RHS becomes nonzero now
    RHSp = new SgNum[ MaxSGLen ];
    VectAssign( RHSp , SgNum( 0 ) , CrrSGLen );
    if( SGBse )
     VectAssignB( RHSp + strt , tmpG1 , SGBse , stp - strt );
    else
     VectAssign( RHSp + strt , tmpG1 , stp - strt );

    RHSxd = HpINF;
    }
   }
  }

 // for each item in the bundle, ask the FiOracle the new information - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < MaxItemN() ; i++ ) {
  if( ! IsThere( i ) )  // for each *existing* item, that is
   continue;

  cIndex_Set SGBse;
  #if( G_IMPLM < 3 )
   SgRow tGi = SubG[ i ];
  #else
   SgRow tGi = tmpG1;
  #endif

  // ask the components [strt , stp) for item i - - - - - - - - - - - - - - -

  Index SGBDim = CrrOrcl->GetGi( tGi + strt , SGBse , i , strt , stp );

  // if they are given in sparse format, densify them - - - - - - - - - - - -

  if( SGBse )
   Densify( tGi , SGBse , SGBDim , stp , strt );

  tGi += strt;

  // add the RHS (if any) - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( RHSp )
   VectSum( tGi , RHSp + strt , stp - strt );

  // copy the information in the row-wise description - - - - - - - - - - - -

  #if( G_IMPLM > 0 )
   #if( HAVE_BTH )
    if( ! ActvVrs )
   #endif
     for( Index j = strt ; j < stp ; )
      TSubG[ j++ ][ i ] = *(tGi++);
  #endif

  }  // end( for( i ) )

 // now have the QP to "digest" the changes - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ChangeQ();  // signal that Q has changed and it has to be recomputed
 #if( HV_NNVAR )
  RecomputeRealAlfa();  // also RealAlfa has to be recomputed
 #endif
 UpdateB();  // because of that, the current factorization and all the rest
             // of the related data structures need be recomputed

 }  // end( QPPenaltyMP::ChgSubG )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

QPPenaltyMP::~QPPenaltyMP()
{
 if( MaxBSize )
  MemDealloc();

 }  // end( ~QPPenaltyMP )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )
 #define SprsV BxdVars || ( ActvVrs && ( ActvVDim < CrrSGLen ) )
#else
 #define SprsV ( ActvVrs && ( ActvVDim < CrrSGLen ) )
 #define MBase2 ActvVrs
#endif

#if( LAZY_Q )

 QuNum QPPenaltyMP::GiTGj( cIndex i , cIndex j )
 {
  if( SprsV )
   if( i == j )
    return( Norm( SubG[ i ] , MBase2 ) );
   else
    return( ScalarProductBB( SubG[ i ] , SubG[ j ] , MBase2 ) );
  else
   if( i == j )
    return( Norm( SubG[ i ] , CrrSGLen ) );
   else
    return( ScalarProduct( SubG[ i ] , SubG[ j ] , CrrSGLen ) );

  }  // end( QPPenaltyMP::GiTGj )

#else  /*-------------------------------------------------------------------*/

 void QPPenaltyMP::GiTG( cIndex i , QuRow Qi , cIndex iMax )
 {
  #if( WHAT_4_Q )  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

   VectAssign( Qi , QuNum( 0 ) , iMax );

   if( SprsV ) {
    cIndex_Set tMB2 = MBase2;
    for( Index h ; ( h = *(tMB2++) ) < InINF ; )
     VectSum( Qi , TSubG[ h ] , TSubG[ h ][ i ] , iMax );
    }
   else {
    SgMat tTSG = TSubG;
    for( Index h = CrrSGLen ; h-- ; tTSG++ )
     VectSum( Qi , *tTSG , (*tTSG)[ i ] , iMax );
    }

  #else  /* WHAT_4_Q == 0- - - - - - - - - - - - - - - - - - - - - - - - - -*/

   SgMat tSG = SubG;
   SgRow SGi = tSG[ i ];
   Index n = 0;

   if( SprsV ) {
    for( ; n < i ; n++ , Qi++ , tSG++ )
     if( IsThere( n ) )
      *Qi = ScalarProductBB( SGi , *tSG , MBase2 );

    *(Qi++) = Norm( *(tSG++) , MBase2 );

    for( ; ++n < iMax ; Qi++ , tSG++ )
     if( IsThere( n ) )
      *Qi = ScalarProductBB( SGi , *tSG , MBase2 );
    }
   else {
    for( ; n < i ; n++ , Qi++ , tSG++ )
     if( IsThere( n ) )
      *Qi = ScalarProduct( SGi , *tSG , CrrSGLen );

    *(Qi++) = Norm( *(tSG++) , CrrSGLen );

    for( ; ++n < iMax ; Qi++ , tSG++ )
     if( IsThere( n ) )
      *Qi = ScalarProduct( SGi , *tSG , CrrSGLen );
    }

  #endif   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  }  // end( QPPenaltyMP::GiTG )

#endif

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

 cSgRow QPPenaltyMP::GiTilde( cIndex i )
 {
  #if( G_IMPLM > 0 )
   return( TSubG[ i ] );
  #else
   return( TSGi( i ) );
  #endif

  }  // end( QPPenaltyMP::GiTilde )

#endif

/*--------------------------------------------------------------------------*/

LMNum QPPenaltyMP::CalculateZ( cIndex h )
{
 #if( WHAT_4_D )
  return( ScalarProduct( Mult , TSubG[ h ] , Base ) );
 #else
  LMNum zh = 0;
  HpRow Mlt = Mult;
  Index_Set Bse = Base;
  for( Index i ; ( i = *(Bse++) ) < InINF ; )
   zh += *(Mlt++) * SubG[ i ][ h ];

  return( zh );
 #endif

 }  // end( QPPenaltyMP::CalculateZ( < one > ) )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::CalculateZ( cIndex_Set Wh , LMRow z )
{
 #if( WHAT_4_D )

  for( Index h ; ( h = *(Wh++) ) < InINF ; )
   z[ h ] = ScalarProduct( Mult , TSubG[ h ] , Base );

 #else

  if( BDim ) {
   cHpRow tM = Mult;
   cIndex_Set tB = Base;
   VectAssignBB( z , SubG[ *tB ] , *tM , Wh );

   for( Index h ; ( h = *(++tB) ) < InINF ; )
    VectSumBB( z , SubG[ h ] , *(++tM) , Wh );
   }
  else
   VectAssignB( z , LMNum( 0 ) , Wh );

 #endif

 }  // end( QPPenaltyMP::CalculateZ( < a set > ) )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::CalculateZ( LMRow z )
{
 if( ActvVrs && ( ActvVDim < CrrSGLen ) )
  QPPenaltyMP::CalculateZ( ActvVrs , z );
 else
  FullZ( z , Base , Mult );

 }  // end( QPPenaltyMP::CalculateZ( < all > ) )

/*--------------------------------------------------------------------------*/

#if( HV_NNVAR )

LMNum QPPenaltyMP::GiTLB( cIndex h , cLMRow l , cIndex_Set lB , cIndex lBd )
{
 #if( G_IMPLM >= 3 )
  if( h != Insrtd ) {
   LMNum res = 0;
   for( Index i ; ( i = *(lB++) ) < InINF ; )
    res += TSubG[ i ][ h ] * l[ i ];

   return( res );
   }
  else
   return( ScalarProductBB( l , tmpG1 , lB ) );
 #else
  return( ScalarProductBB( l , SubG[ h ] , lB ) );
 #endif
 }

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::GiTLB( HpRow gtlb , cLMRow l , cIndex_Set lB ,
			 cIndex lBd , const bool add )
{
 #if( WHAT_4_Q )  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  for( Index i ; ( i = *(lB++) ) < InINF ; )
   if( ABS( l[ i ] ) > Eps< HpNum >() )
    VectSum( gtlb , TSubG[ i ] , add ? l[ i ] : - l[ i ] , MaxItemN() );

 #else  /* WHAT_4_Q == 0- - - - - - - - - - - - - - - - - - - - - - - - - -*/

  if( add ) {
   for( Index n = 0 ; n < MaxItemN() ; n++ , gtlb++ )
    if( IsThere( n ) )
     *gtlb += ScalarProductBB( l , SubG[ n ] , lB );
   }
  else
   for( Index n = 0 ; n < MaxItemN() ; n++ , gtlb++ )
    if( IsThere( n ) )
     *gtlb -= ScalarProductBB( l , SubG[ n ] , lB );

 #endif   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 }

#endif

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void QPPenaltyMP::ComputeRHSxd( void )
{
 if( RHSxd == HpINF ) {  // still have to compute RHS * d
  if( RHSp )
   #if( HV_NNVAR )
    RHSxd = DPerG( RHSp );
   #else
    {
     CompleteZ( false );

     if( ActvVrs && ( ActvVDim < CrrSGLen ) )
      RHSxd = ScalarProduct( dir , RHSp , ActvVrs );
     else
      RHSxd = ScalarProduct( dir , RHSp , CrrSGLen );
     }
   #endif
  else
   RHSxd = 0;
  }
 }  // end( ComputeRHSxd )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::Computedir( bool Fulld )
{
 if( dirCmptd < ( Fulld ? 2 : 1 ) ) {
  // calculate the entries of z*: they sit in the internal data structure - -
  // of BMinQuad if HV_NNVAR, and in dir[] otherwise- - - - - - - - - - - - -

  CompleteZ( Fulld );

  // now calculate the actual d - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( HV_NNVAR )
   ReadD( dir , Fulld ? CrrSGLen : 0 );
  #else
   if( Fulld )
    VectScale( dir , - tCurr , CrrSGLen );
   else {
    cIndex_Set tAV = ActvVrs;
    for( Index h ; ( h = *(tAV++) ) < InINF ; )
     dir[ h ] *= - tCurr;
    }

   // note that if BMinQuad is not used, z* and d* "share the same memory":
   // hence, ZCmptd > 1 ==> dirCmptd == 0 && dirCmptd > 0 ==> ZCmptd == 1
   ZCmptd = 1;
  #endif

  dirCmptd = Fulld ? 2 : 1;
  }
 }  // end( Computedir )

/*--------------------------------------------------------------------------*/

Index QPPenaltyMP::CheckBCopy( void )
{
 // check if the item is identical to any other in the bundle - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index IsIde = InINF;

 if( ChkIde )
  for( Index i = 0 ; i < MaxItemN() ; i++ ) {
   #if( ADD_CNST )
    if( G1IsSubg ) {
     if( ! IsASubG( i ) )
      continue;
     }
    else
     if( ! IsAConst( i ) )
      continue;
   #else
    if( ! IsThere( i ) )
     continue;
   #endif

   #if( G_IMPLM < 3 )
    if( EqualVect( tmpG1 , SubG[ i ] , CrrSGLen ) ) {
     IsIde = i;
     break;
     }
   #else
    SgRow tGi = tmpG1;
    Index j = CrrSGLen;
    for( SgMat tTSG = TSubG ; j ; j-- )
     if( (*(tTSG++))[ i ] != *(tGi++) )
      break;

    if( ! j ) {
     IsIde = i;
     break;
     }
   #endif
   }

 return( IsIde );

 }  // end( CheckBCopy )

/*--------------------------------------------------------------------------*/

#if( HAVE_BTH )

void QPPenaltyMP::GetNewGTRow( cIndex i )
{
 if( ! BFCntr )
  for( Index j = BFCntr = std::min( MaxSGLen - NrAllRws , GTSltSize ) ; j-- ; )
   Buffer[ j ] = new SgNum[ MaxBSize ];

 TSubG[ i ] = Buffer[ --BFCntr ];

 NrAllRws++;
 }

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::SetGTRow( cIndex i )
{
 // copy the entries from the columns into TSubG[ i ]

 SgMat sgt = SubG;     
 SgRow tTSGk = TSubG[ i ];
 for( Index j = MaxItemN() ; j-- ; )
  *(tTSGk++) = (*(sgt++))[ i ];
 }

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::PutGTRow( cIndex i )
{
 Buffer[ BFCntr++ ] = TSubG[ i ];
 TSubG[ i ] = 0;
 NrAllRws--;
 }

#endif

/*--------------------------------------------------------------------------*/

#if( G_IMPLM < 1 )

SgRow QPPenaltyMP::TSGi( cIndex i )
{
 SgMat sgt = SubG;
 SgRow tTSGk = TSGk;
 for( Index j = MaxItemN() ; j-- ; )
  *(tTSGk++) = (*(sgt++))[ i ];

 return( TSGk );
 }

#endif

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::FullZ( LMRow z , cIndex_Set tB , cHpRow tM )
{
 // computes z for the set of multipliers given by tB and tM, comprised the
 // entries corresponding to variables that are not defined in the QP
 // subproblem; if HAVE_BTH the row-wise representation is preferred, but it
 // can be used only if no L.V.G. is being used, since in this case not all
 // the entries are in principle available row-wise

 #if( WHAT_4_D && HAVE_BTH )   /*- - - - - - - - - - - - - - - - - - - - -*/
  if( ! ActvVrs )
 #endif
  #if( WHAT_4_D )
  {
   SgMat tTSG = TSubG + CrrSGLen;
   for( z += CrrSGLen ; tTSG > TSubG ; )
    *(--z) = ScalarProduct( tM , *(--tTSG) , tB );
   }
  #endif
 #if( WHAT_4_D && HAVE_BTH )
  else
 #endif
  #if( ( ! WHAT_4_D ) || HAVE_BTH )  /*- - - - - - - - - - - - - - - - - -*/
   if( BDim ) {
    VectAssign( z , SubG[ *tB ] , *tM , CrrSGLen );

    for( Index h ; ( h = *(++tB) ) < InINF ; )
     VectSum( z , SubG[ h ] , *(++tM) , CrrSGLen );
    }
   else
    VectAssign( z , LMNum( 0 ) , CrrSGLen );

  #endif   /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 }  // end( FullZ )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::CompleteZ( bool Fullz )
{
 // if the active set is either not defined or full, Fullz == false has no
 // sense
 if( ( ! ActvVrs ) || ( ActvVDim == CrrSGLen ) )
  Fullz = true;

 #if( ! HV_NNVAR )
  // if dir[] contains d*, revert it to z*- - - - - - - - - - - - - - - - - -
  // note that if BMinQuad is not used, z* and d* "share the same memory":
  // hence, ZCmptd > 1 ==> dirCmptd == 0 && dirCmptd > 0 ==> ZCmptd == 1

  LMRow zz = dir;

  if( dirCmptd == 2 ) {
   VectScale( zz , 1 / ( - tCurr ) , CrrSGLen );
   dirCmptd = 0;
   ZCmptd = 4;
   }
  else
   if( dirCmptd == 1 ) {
    VectScale( zz , 1 / ( - tCurr ) , ActvVrs );
    dirCmptd = 0;
    ZCmptd = 3;
    }
 #else
  LMRow zz = SetD();
 #endif

 if( ZCmptd >= ( Fullz ? 4 : 3 ) )
  return;

 if( ZCmptd < 1 )
  throw( NDOException( "QPPenaltyMP::MakeZ: QP not (correctly) solved" ) );

 #if( HAVE_BTH )
  // if HAVE_BTH and some variables are not active, TSubG[ h ]- - - - - - - -
  // is not defined for the inactive h, hence all z* must be computed

  if( Fullz && ActvVrs && ( ActvVDim < CrrSGLen ) ) {
   FullZ( zz , MinQuad::ReadBase() , MinQuad::ReadMult() );
   ZCmptd = 4;
   #if( ! HV_NNVAR )
    dirCmptd = 0;
   #endif
   return;
   }
 #endif

 #if( WHAT_4_D )
  // the row-wise form allows partial recalculations- - - - - - - - - - - - -

  if( ZCmptd > 1 ) {
   cHpRow tM = MinQuad::ReadMult();
   cIndex_Set tB = MinQuad::ReadBase();

   #if( HV_NNVAR && ( LAZY_D == 2 ) )
    if( ZCmptd == 2 )  // - - - - - - - - - - - - - - - - - - - - - - - - - -
     if( Fullz ) {
      // compute the entries of z corresponding to variables that are either
      // undefined or UC

      const char NNP = NNVar() | IsVar();
      for( Index h = 0 ; h < CrrSGLen ; h++ )
       if( ( GetVar( h ) & NNP ) != NNP )
        SetD( h , ScalarProduct( tM , TSubG[ h ] , tB ) );

      ZCmptd = 4;
      }
     else {
      // compute the entries of z corresponding to defined UC variables

      cIndex_Set tIAV = InActiveVars();
      for( Index h ; ( h = *(tIAV++) ) < InINF ; )
       if( ! ( GetVar( h ) & NNVar() ) )
        SetD( h , ScalarProduct( tM , TSubG[ h ] , tB ) );

      ZCmptd = 3;
      }
    else  // ZCmptd == 3 - - - - - - - - - - - - - - - - - - - - - - - - - -
   #endif
    {
     cIndex_Set tAV = ActvVrs;
     for( Index h = 0 ; h < CrrSGLen ; h++ )
      if( *tAV == h )
       tAV++;
      else
       zz[ h ] = ScalarProduct( tM , TSubG[ h ] , tB );

     ZCmptd = 4;
     }
   }
  else  // ZCmptd == 1
 #endif
  {
   // if all z* must be computed or only the column-wise form - - - - - - - -
   // is available, nothing terribly smart can be done

   if( Fullz ) {
    FullZ( zz , MinQuad::ReadBase() , MinQuad::ReadMult() );
    ZCmptd = 4;
    }
   else {
    QPPenaltyMP::CalculateZ( ActvVrs , zz );
    ZCmptd = 3;
    }
   }

 #if( ! HV_NNVAR )
  dirCmptd = 0;
 #endif

 }  // end( CompleteZ )

/*--------------------------------------------------------------------------*/

#if( ADD_CNST )

void QPPenaltyMP::ComputeZBase( void )
{
 ZBDim = 0;
 cHpRow tM = MinQuad::ReadMult();
 cIndex_Set tB = MinQuad::ReadBase();
 for( Index i ; ( i = *(tB++) ) < InINF ; ++tM )
  if( ! IsAConst( i ) ) {
   ZBase[ ZBDim ] = i;
   ZMult[ ZBDim++ ] = *tM;
   }
 }  // end( ComputeZBase )

#endif

/*--------------------------------------------------------------------------*/

#define MPLOG( x ) if( MPLog ) *MPLog << x;

void QPPenaltyMP::ChgRHS( cIndex strt , cIndex stp )
{
 if( ! MaxItemN() ) {  // special case: the bundle is empty - - - - - - - - -
                       // - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( strt == 0 ) && ( stp == CrrSGLen ) ) {
   // even more special case: changing all of the RHS - - - - - - - - - - - -
   if( G1Perz ) {          // ... which is not all-0
    Swap( tmpG1 , RHSp );  // put it in place

    if( ! tmpG1 )          // if necessary, reallocate tmpG1
     tmpG1 = new SgNum[ MaxSGLen ];
    }
   else {                  // the RHS is all-0
    delete[] RHSp;
    RHSp = nullptr;
    }
   }
  else  // changing only a part of the RHS (in an empty bundle) - - - - - - -
   if( RHSp )              // there was a (nonzero) RHS set already
    if( G1Perz )           // and there are more nonzeroes
     VectAssign( RHSp + strt , tmpG1 , stp - strt );
    else                   // all the new part is zero
     VectAssign( RHSp + strt , SgNum( 0 ) , stp - strt );
   else                    // the previous RHS was all-zero
    if( G1Perz ) {         // but this is not
     RHSp = new SgNum[ MaxSGLen ];                // first allocate it
     VectAssign( RHSp , SgNum( 0 ) , CrrSGLen );  // and zero it all
     VectAssign( RHSp + strt , tmpG1 , stp - strt );
     }
    // else the RHS was all-zero and a part has set to all-zero, pretty
    // little to do ...

  RHSxd = HpINF;
  return;                  // nothing else to do
  }

 // compute the deltaRHS = new RHS - old RHS in tmpG1 - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // meanwhile, ensure that the new RHS is saved into RHSp

 RHSxd = HpINF;  // we assume that something will change ...

 if( RHSp )       // if a RHS is set already- - - - - - - - - - - - - - - - -
  if( G1Perz )    // and the new RHS is not all-0
   if( ( strt == 0 ) && ( stp == CrrSGLen ) ) {  // changing all the RHS
    Swap( tmpG1 , RHSp );  // now RHSp = new RHS, tmpG1 = old RHS
    for( Index i = 0 ; i < CrrSGLen ; i++ )
     tmpG1[ i ] = RHSp[ i ] - tmpG1[ i ];
     }
   else                                          // changing only a part
    for( Index i = 0 ; i < stp - strt ; i++ ) {
     cSgNum nRHSi = tmpG1[ i ];
     tmpG1[ i ] -= RHSp[ strt + i ];
     RHSp[ strt + i ] = nRHSi;
     }
  else            // but the new RHS is all-0 ==> deltaRHS = - oldRHS
   if( ( strt == 0 ) && ( stp == CrrSGLen ) ) {  // changing all the RHS
    for( Index i = 0 ; i < CrrSGLen ; i++ )
     tmpG1[ i ] = - RHSp[ i ];

    delete[] RHSp;  // deallocate RHSp
    RHSp = 0;       // the RHS is all-0
    }
   else                                          // changing only a part
    for( Index i = 0 ; i < stp - strt ; i++ ) {
     tmpG1[ i ] = - RHSp[ strt + i ];
     RHSp[ strt + i ] = 0;
     }
 else             // the RHS is currently all-0 - - - - - - - - - - - - - - -
                  // which implies deltaRHS = newRHS - 0 = tmpG1 already
  if( G1Perz ) {  // and the new RHS is not all-0
   RHSp = new SgNum[ MaxSGLen ];  // first allocate the RHSp
   if( ( strt == 0 ) && ( stp == CrrSGLen ) )  // changing all the RHS
    VectAssign( RHSp , tmpG1 , CrrSGLen );     // RHS = newRHS
   else {                                      // changing only a part
    HpRow tRHS = RHSp;
    for( Index i = strt ; i-- ; )
     *(tRHS++) = 0;
    for( Index i = 0 ; i < stp - strt ; i++ )
     *(tRHS++) = tmpG1[ i ];
    for( Index i = CrrSGLen - stp ; i-- ; )
     *(tRHS++) = 0;
    }
   }
  else {          // but the new RHS is all-0
   MPLOG( "QPPenaltyMP::ChgRHS: RHS changed from 0 to 0?" );
   return;        // well, not much to do, is there?
   }

 // now update all items in the bundle- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( G_IMPLM < 3 )
  // update the column-wise representation- - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < MaxItemN() ; i++ )
   if( IsThere( i ) )
    VectSum( SubG[ i ] + strt , tmpG1 , stp - strt );
 #endif

 #if( G_IMPLM > 0 )
  // update the row-wise representation - - - - - - - - - - - - - - - - - - -

  #if( HAVE_BTH )
   if( ActvVrs && ( ActvVDim < CrrSGLen ) ) {
    cIndex_Set AVt = ActvVrs;
    while( *AVt < strt )
     AVt++;
    for( Index h ; ( h = *(AVt++) ) < stp ; )
     VectSum( TSubG[ h ] , tmpG1[ h - strt ] , MaxItemN() );
    }
   else
  #endif
    for( Index i = 0 ; i < stp - strt ; i++ )
      VectSum( TSubG[ strt + i ] , tmpG1[ i ] , MaxItemN() );
 #endif

 // now have the QP to "digest" the changes - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 /* Note: there surely are ways to update Q & L more effectively than how
    it is done now, but there will (perhaps) be time for that. */

 ChangeQ();  // signal that Q has changed and it has to be recomputed
 #if( HV_NNVAR )
  RecomputeRealAlfa();  // also RealAlfa has to be recomputed
 #endif
 UpdateB();  // because of that, the current factorization and all the rest
             // of the related data structures need be recomputed

 }  // end( ChgRHS )

/*--------------------------------------------------------------------------*/

#define AREDIFF( x , y ) ABS( x - ( y ) ) > \
                               EpsilonD() * std::max( ABS( x ) , HpNum( 1 ) )
#if( HV_NNVAR )

void QPPenaltyMP::CheckBounds( void )
{
 #if( CHECK_DS & 1 )
  cLMRow lb = LowerBounds();
  #if( TWOSIDED )
   cLMRow ub = UpperBounds();
  #endif

  // access the Bundle solver through the FiOracle

  NDOSolver *NDOS = CrrOrcl->GetNDOSolver();

  // read the current point in the Bundle solver

  Index CPBd;
  cIndex_Set CPB;
  cLMRow CP = NDOS->ReadSol( CPB , CPBd );

  // now check the bounds

  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   LMNum CPi;
   if( CPB )
    if( *CPB == i ) {
     CPB++;
     CPi = *(CP++);
     }
    else
     CPi = 0;
   else
    CPi = *(CP++);

   if( lb[ i ] > - LMINF ) {
    if( CrrOrcl->GetUC( i ) )
     MPLOG( endl << "Error: " << i << " is not constrained" );
    else
     if( AREDIFF( CPi , - lb[ i ] ) )
      MPLOG( endl << "Error: lb[ " << i << " ] = " << lb[ i ]
	          << " != - " << CPi );
    }
   else
    if( ! CrrOrcl->GetUC( i ) )
     MPLOG( endl << "Error: " << i << " is constrained" );

   #if( TWOSIDED )
    if( ub[ i ] < LMINF ) {
     if( CrrOrcl->GetUB( i ) == LMINF )
      MPLOG( endl << "Error: " << i << " is not constrained" );
     else
      if( AREDIFF( ub[ i ] , CrrOrcl->GetUB( i ) - CPi ) )
       MPLOG( endl << "Error: ub[ " << i << " ] = " << ub[ i ]
	           << " != - " << CrrOrcl->GetUB( i ) - CPi );
     }
    else
     if( CrrOrcl->GetUB( i ) < LMINF )
      MPLOG( endl << "Error: " << i << " is constrained" );
   #endif

   }  // end( for( i ) )
 #endif

 }  // end( CheckBounds )

#endif

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::CheckSG( void )
{
 #if( CHECK_DS & 2 )
  // check the RHS- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // get the RHS again from the oracle (store it in tmpG1[])- - - - - - - - -

  cIndex_Set SGBse;
  Index SGBDm = CrrOrcl->GetGi( tmpG1 , SGBse , MaxBSize );

  if( SGBse )
   Densify( tmpG1 , SGBse , SGBDm , CrrSGLen );  // in "dense" format

  if( RHSp ) {  // if it is nonzero - - - - - - - - - - - - - - - - - - - - -
   for( Index j = 0 ; j < CrrSGLen ; j++ )
    if( AREDIFF( tmpG1[ j ] , RHSp[ j ] ) )
     MPLOG( endl << "Error in RHS[ " << j << " ] = "
	         << tmpG1[ j ] - RHSp[ j ] );
   }
  else        // if it is all-zero- - - - - - - - - - - - - - - - - - - - - -
   for( Index j = 0 ; j < CrrSGLen ; j++ )
    if( ABS( tmpG1[ j ] ) > EpsilonD() )
     MPLOG( endl << "Error in RHS[ " << j << " ] = " << tmpG1[ j ] );

  // check SubG[] / TSubG[] - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < MaxItemN() ; i++ )
   if( IsThere( i ) ) {  // item `i' is kept in the bundle
    // get G[ i ] again from the oracle (store it in tmpG1[]) - - - - - - - -

    SGBDm = CrrOrcl->GetGi( tmpG1 , SGBse , i );

    if( SGBse )
     Densify( tmpG1 , SGBse , SGBDm , CrrSGLen );  // in "dense" format

    if( RHSp )
     VectSum( tmpG1 , RHSp , CrrSGLen );  // ... and including the RHS

    #if( G_IMPLM >= 3 )
     Insrtd = i;  // keep track of the contents of tmpG1
    #else
     // check SubG[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     for( Index j = 0 ; j < CrrSGLen ; j++ )
      if( AREDIFF( tmpG1[ j ] , SubG[ i ][ j ] ) )
       MPLOG( endl << "Error in SubG[ " << i << " ][ " << j << " ] = "
	           << tmpG1[ j ] - SubG[ i ][ j ] );
    #endif

    #if( G_IMPLM > 0 )
     // check TSubG[] - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     #if( HAVE_BOTH )
      if( ActvVrs && ( ActvVDim < CrrSGLen ) ) {
       cIndex_Set tAV = ActvVrs;
       for( Index j ; ( j = *(tAV++) ) < InINF ; )
	if( AREDIFF( tmpG1[ j ] , TSubG[ j ][ i ] ) )
	 MPLOG( endl << "Error in TSubG[ " << j << " ][ " << i << " ] = "
		     << tmpG1[ j ] - TSubG[ j ][ i ] );
       }
      else
     #else
       for( Index j = 0 ; j < CrrSGLen ; j++ )
	if( AREDIFF( tmpG1[ j ] , TSubG[ j ][ i ] ) )
	 MPLOG( endl << "Error in TSubG[ " << j << " ][ " << i << " ] = "
		     << tmpG1[ j ] - TSubG[ j ][ i ] );
     #endif
    #endif
    }

 #endif

 }  // end( CheckSG )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::CheckGT( void )
{
 #if( CHECK_DS & 4 )
  #if( HV_NNVAR )
   cLMRow tmpdir = Readd( true );  // read the current dir
  #else
   // compute z* anew

   LMRow tmpdir = new LMNum[ CrrSGLen ];
   FullZ( tmpdir , MinQuad::ReadBase() , MinQuad::ReadMult() );

  // check against existing z*, if any

   if( ZCmptd == 4 )
    for( Index j = 0 ; j < CrrSGLen ; j++ ) {
     if( AREDIFF( tmpdir[ j ] , dir[ j ] ) )
      MPLOG( endl << "Error in z*[ " << j << " ] = "
	          << tmpdir[ j ] - dir[ j ] );
     }
   else
    if( ZCmptd == 3 ) {
     cIndex_Set tAV = ActvVrs;
     for( Index j ; ( j = *(tAV++) ) < InINF ; )
      if( AREDIFF( tmpdir[ j ] , dir[ j ] ) )
       MPLOG( endl << "Error in z*[ " << j << " ] = "
	           << tmpdir[ j ] - dir[ j ] );
     }
  #endif

  // now check scalar products

  for( Index i = 0 ; i < MaxItemN() ; i++ )
   if( IsThere( i ) ) {
    #if( G_IMPLM < 3 )
     SgRow SGi = SubG[ i ];
    #else
     for( Index j = 0 ; j < CrrSGLen ; j++ )
      tmpG1[ j ] = TSubG[ j ][ i ];

     SgRow SGi = tmpG1;
     Insrtd = i;
    #endif

    HpNum GTdi;
    if( ActvVrs && ( ActvVDim < CrrSGLen ) )
     GTdi = ScalarProductBB( tmpdir , SGi , ActvVrs );
    else
     GTdi = ScalarProduct( tmpdir , SGi , CrrSGLen );

    #if( HV_NNVAR )
     // the scalar product is with d* = - t z*, but ReadGTz( i ) returns
     // the opposite value, so scale it accordingly (and change its sign)
     GTdi *= - 1 / tCurr;
    #endif

    if( AREDIFF( GTdi , ReadGTz( i ) ) )
     MPLOG( endl << "Error in G[ " << i << " ] * d = "
	         << GTdi - ReadGTz( i ) );
    }

  #if( ! HV_NNVAR )
   // turn z* into d* and check against the existing one, if any 

   if( dirCmptd == 1 ) {
    VectScale( tmpdir , - tCurr , CrrSGLen );

    if( dirCmptd == 2 )
     for( Index j = 0 ; j < CrrSGLen ; j++ ) {
      if( AREDIFF( tmpdir[ j ] , dir[ j ] ) )
       MPLOG( endl << "Error in d*[ " << j << " ] = "
	           << tmpdir[ j ] - dir[ j ] );
      }
    else {
     cIndex_Set tAV = ActvVrs;
     for( Index j ; ( j = *(tAV++) ) < InINF ; )
      if( AREDIFF( tmpdir[ j ] , dir[ j ] ) )
       MPLOG( endl << "Error in d*[ " << j << " ] = "
	           << tmpdir[ j ] - dir[ j ] );
     }
    }

   // cleanup temporary data structure
   delete[] tmpdir;
  #endif
 #endif

 }  // end( CheckGT )

/*--------------------------------------------------------------------------*/

void QPPenaltyMP::MemDealloc( void )
{
 delete[] RHSp;

 #if( ADD_CNST )
  delete[] ZMult;
  delete[] ZBase;
 #endif

 #if( G_IMPLM < 3 )
  for( ; DimMinQuad-- ; )
   delete[] SubG[ DimMinQuad ];

  delete[] SubG;
 #endif

 delete[] tmpG1;

 #if( G_IMPLM > 0 )
  #if( G_IMPLM < 3 )
   for( ; BFCntr-- ; )
    delete[] Buffer[ BFCntr ];

   delete[] Buffer;
  #endif

  if( TSubG ) {
   for( Index i = MaxSGLen ; i-- ; )
    delete[] TSubG[ i ];

   delete[] TSubG;
   }
 #endif

 #if( G_IMPLM < 1 )
  delete[] TSGk;
 #endif

 delete[] dir;

 }  // end( MemDealloc )

/*--------------------------------------------------------------------------*/
/*-------------------------- End File QPPnltMP.C ---------------------------*/
/*--------------------------------------------------------------------------*/
