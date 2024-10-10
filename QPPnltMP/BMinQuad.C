/*--------------------------------------------------------------------------*/
/*--------------------------- File BMinQuad.C ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                                                                      --*/
/*-- Implementation of the BTT algorithm for solving the box-constrained  --*/
/*-- Quadratic Problems arising as descent direction finding subproblems  --*/
/*-- within (box)Constrained Bundle algorithms. Uses as a subroutine the  --*/
/*-- TT algorithm for linearly constrained problems without bounds.       --*/
/*--                                                                      --*/
/*--                         Antonio Frangioni                            --*/
/*--                    Dipartimento di Informatica                       --*/
/*--                        Universita' di Pisa                           --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BMinQuad.h"

#include "OPTvect.h"

#include <assert.h>

#include <algorithm>

#include <cmath>

/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--                                                                      --*/
/*--      Some small macro definitions, used throughout the code.         --*/
/*--                                                                      --*/
/*--------------------------------------------------------------------------*/

#define CHECK_DS 0

/* If CHECK_DS > 0, various data structures are checked for correctness
   during the run of the algorithm, tremendously slowing down the algorithm
   but allowing to debug the thing.

   What data structures are checked is coded bit-wise in CHECK_DS:

    bit 0 (+ 1)  =>  Base2[], MBase2[] and GS[] are checked
    bit 1 (+ 2)  =>  RealAlfa[] is checked

   CHECK_DS > 0 forces asserts() do work within this unit. */

#if( CHECK_DS )
 #ifdef NDEBUG
  #undef NDEBUG
 #endif
#endif

#if( TWOSIDED )
 #define lb( h ) lb[ h ]
#else
 #define lb( h ) bounds[ h ]
#endif

#define size( x ) eDir * std::max( ABS( x ) , LMNum( 1 ) )

#if( ! SIGNAL_B2CHG )
 #define B2HasChgd()
#endif

#if( ! SIGNAL_MBCHG )
 #define MBHasChgd()
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace MinQuad_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

static cHpNum MaxEpsD = 1e-3;   // Maximum, ...
static cHpNum IntEpsD = 1e-6;   // and Initial value for eD

static cHpNum MFactor = 10;     // eD is increased by multiplying
                                // it by this factor

/* Meaning of the bits in the char describing a variable:

   bit 0 ==> 0 if it is UC, 1 if it is NN;
   bit 1 ==> 0 if it is non-existing, 1 if it is existing;
   bit 2 ==> 0 if it is in MBase2, 1 if it is in Base2;
   bit 3 ==> 0 if it is in Base2 at the LB, 1 if it is in Base2 at the UB.
   bit 4 ==> 0 if it is "normal", 1 if it is "taboo". */

static const char kNNN   =  1;  // a Non-existing NN variable
static const char kEUC   =  2;  // an Existing UC variable (in MBase2)
static const char kENN   =  3;  // an Existing NN variable in MBase2
static const char kLBC   =  7;  // an Existing NN variable in Base2 at its LB
#if( TWOSIDED )
 static const char kUBC  = 15;  // an Existing NN variable in Base2 at its UB
#endif

static const char kIsNN  =  1;  // ( i &  1 ) <=> i is a NN variable
static const char kIsIn  =  2;  // ( i &  2 ) <=> i is existing
static const char kInB2  =  4;  // ( i &  4 ) <=> i is in Base2
#if( TWOSIDED )
 static const char kIsUB =  8;  // ( i &  8 ) <=> i is in Base2 at the UB
#endif
static const char kIsTb  = 16;  // ( i & 16 ) <=> i is "taboo"

static cIndex InINF = Inf< Index >();
static cHpNum HpINF = Inf< HpNum >();
static cHpNum HpEps = Eps< HpNum >();

/*--------------------------------------------------------------------------*/
/*-------------------- IMPLEMENTATION OF BMinQuad  -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTRUCTOR -- -----------------------------*/
/*--------------------------------------------------------------------------*/

BMinQuad::BMinQuad( void )
          :
          MinQuad()
{
 MaxVarAdd = MaxVarRmv = InINF;

 #if( LOG_BMQ )
  BMQLog = &clog;

  #if( LOG_BMQ > 1 ) 
   BCalls = BSccss = 0;
   SumAverages = 0;
  #endif
 #endif

 BMQt = 0;

 }  // end( BMinQuad )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BMinQuad::SetMaxDim( Index m , Index n , Index SDim )
{
 if( MaxBDim )   // this is not the first call: - - - - - - - - - - - - - - -
  MemDealloc();  // deallocate everything - - - - - - - - - - - - - - - - - -

 MinQuad::SetMaxDim( m , n , SDim );  // "throw" the method of the base class

 if( MaxBDim ) {    // m != 0: allocate everything- - - - - - - - - - - - - -
  SpaceDim = SDim;  //- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  RealAlfa = new HpNum[ m ];

  GS       = new char[ SpaceDim ];
  Base2    = new Index[ SpaceDim + 2 ];
  MBase2   = Base2 + SpaceDim + 1;
  MvdVars  = new Index[ SpaceDim + 2 ];

  di       = new LMNum[ SpaceDim ];
  #if( ! CNDVD_TMP )
   tmpdi   = new LMNum[ SpaceDim ];
  #endif
  bounds   = new LMNum[ SpaceDim ];

  #if( TWOSIDED )
   lb      = new LMNum[ SpaceDim ];
   ub      = new LMNum[ SpaceDim ];
  #endif

  for( Index i = SpaceDim ; i-- ; ) {
   GS[ i ] = kNNN;           // by default, all variables are NN
   lb( i ) = 0;              // with 0 LB
   #if( TWOSIDED )
    ub[ i ] = Inf< LMNum >();  // and +INF UB
   #endif
   }

  *MBase2 = *Base2 = InINF;

  }  // end( if( m != 0 ) )- - - - - - - - - - - - - - - - - - - - - - - - - -

 // variables initialization - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 NNStop = TLDim = MB2Dim = B2Dim = 0;
 eD = IntEpsD;
 Bf = HpINF;

 #if( ! BEXACT )
  bNorm = 0;
 #endif

 }  // end( BMinQuad::SetMaxDim )

/*--------------------------------------------------------------------------*/

void BMinQuad::SetEpsilonD( HpNum NeweD )
{
 eD = NeweD ? NeweD : IntEpsD;
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR ADDING / REMOVING / CHANGING DATA -------------*/
/*--------------------------------------------------------------------------*/

void BMinQuad::AddSubGrad( cIndex n , cHpNum alfan )
{
 Bf = HpINF;  // the objective function is allowed to increase

 RealAlfa[ n ] = alfan;
 cHpNum GnTLB = GiTLB( n , bounds , Base2 , B2Dim );

 MinQuad::AddSubGrad( n , alfan - GnTLB );

 MinQuad::SetGTz( n , GnTLB / PrvsTi );

 }  // end( BMinQuad::AddSubGrad )

/*--------------------------------------------------------------------------*/

void BMinQuad::AddConstr( cIndex n , cHpNum alfan )
{
 Bf = HpINF;  // the objective function is allowed to increase

 RealAlfa[ n ] = alfan;
 cHpNum GnTLB = GiTLB( n , bounds , Base2 , B2Dim );

 MinQuad::AddConstr( n , alfan - GnTLB );

 MinQuad::SetGTz( n , GnTLB / PrvsTi );

 }  // end( BMinQuad::AddConstr )

/*--------------------------------------------------------------------------*/

void BMinQuad::ChangeAlfa( cHpNum DeltaAlfa )
{
 for( Index i = 0 ; i < NxtBIdx ; i++ )
  if( IsASubG( i ) )
   RealAlfa[ i ] += DeltaAlfa;

 MinQuad::ChangeAlfa( DeltaAlfa );

 Bf = HpINF;  // the objective function is allowed to increase

 }  // end( ChangeAlfa( HpNum , HpRow ) )

/*--------------------------------------------------------------------------*/

void BMinQuad::MoveAlongD( cHpNum Tau , cHpNum DeltaFi )
{
 // take care of "active" constraints - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index h;
 cIndex_Set o = Base2;

 if( Tau >= PrvsTi ) {  //- - - - - - - - - - - - - - - - - - - - - - - - - -
  // if Tau >= PrvsTi, the components of the direction which hit the bounds
  // for PrvsTi (precisely the ones in Base2[]) still hit the bounds for Tau

  for( ; ( h = *(o++) ) < InINF ; ) {
   #if( TWOSIDED )
   if( lb[ h ] > - Inf< LMNum >() )
     lb[ h ] -= bounds[ h ];

   if( ub[ h ] < Inf< LMNum >() )
     ub[ h ] -= bounds[ h ];
   #endif

   bounds[ h ] = 0;
   }

  #if( ! BEXACT )
   bNorm = 0;
  #endif
  }
 else {  // Tau < PrvsTi- - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the components of the direction which hit the bounds for PrvsTi do not
  // do that if a shorter step is taken

  // if Tau == 0, nothing really changes
  if( Tau > 0 ) {
   cHpNum tauR = Tau / PrvsTi;
   cHpNum tauRM1 = 1 - tauR; 

   for( ; ( h = *(o++) ) < InINF ; ) {
    #if( TWOSIDED )
     if( lb[ h ] > - Inf< LMNum >() )
      lb[ h ] -= tauR * bounds[ h ];

     if( ub[ h ] < Inf< LMNum >() )
      ub[ h ] -= tauR * bounds[ h ];
    #endif

    bounds[ h ] *= tauRM1;
    }

   #if( ! BEXACT )
    bNorm *= ( tauRM1 * tauRM1 );
   #endif
   }
  }  // end else( Tau < PrvsTi )- - - - - - - - - - - - - - - - - - - - - - -

 // take care of "inactive" constraints - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // if Tau == 0, nothing really changes
 if( Tau > 0 ) {
  for( o = MBase2 ; ( h = *(o++) ) < InINF ; ) {
   #if( TWOSIDED )
    cLMNum dh = di[ h ] * Tau;

    if( ub[ h ] < + Inf< LMNum >() )
     ub[ h ] += dh;

    if( lb[ h ] > - Inf< LMNum >() )
     lb[ h ] += dh;
   #else
    if( GS[ h ] & kIsNN )
     bounds[ h ] += di[ h ] * Tau;
   #endif
   }
  }

 // take care of Alfa[] and RealAlfa[]- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // we have that Alfa[] = RealAlfa[] - G_{Xsi} l_{Xsi}, and that both
 // RealAlfa[] and l_{Xsi} change; in particular,
 //  l_{Xsi} := l_{Xsi} * ( 1 - Tau / t )
 //  RealAlfa := RealAlfa[] - ( Tau / t ) * G * d + DeltaFi
 // The result is that MinQuad correctly updates Alfa[], i.e., Alfa[] at
 // the end of MinQuad::MoveAlongD() correctly takes into account the
 // term G_{Xsi} l_{Xsi} for the new l_{Xsi}, while RealAlfa[] has to
 // be kept updated
 //
 // however, if Tau == 0 then G_{Xsi} l_{Xsi} does not change, and since the
 // change of Alfa[] is additive then also that of RealAlfa is

 if( ( Tau > 0 ) && ( Tau != PrvsTi ) ) {
  // RealAlfa := G_{Xsi} l_{Xsi}
  VectSubtract( RealAlfa , Alfa , NxtBIdx );

  // RealAlfa := ( 1 - Tau / t ) * G_{Xsi} l_{Xsi}
  VectScale( RealAlfa , 1 - Tau / PrvsTi , NxtBIdx );
  }

 MinQuad::MoveAlongD( Tau , DeltaFi );

 if( Tau > 0 ) {
  if( Tau == PrvsTi )
   VectAssign( RealAlfa , Alfa , NxtBIdx );
  else
   VectSum( RealAlfa , Alfa , NxtBIdx );
  }

 CheckRA();

 Bf = HpINF;  // the objective function is allowed to increase

 }  // end( BMinQuad::MoveAlongD )

/*--------------------------------------------------------------------------*/

void BMinQuad::ReadAlfa( HpRow NewAlfa )
{
 VectAssign( NewAlfa , RealAlfa , NxtBIdx );

 }  // end( BMinQuad::ReadAlfa( HpRow ) )

/*--------------------------------------------------------------------------*/

void BMinQuad::SetGTz( cIndex i , cHpNum GTzi )
{
 MinQuad::SetGTz( i , GTzi + MinQuad::ReadGTz( i ) );
 }

/*--------------------------------------------------------------------------*/

void BMinQuad::ChangeBounds( void )
{
 // recompute bNorm - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bNorm = 0;
 #if( TWOSIDED )
  Index Moved = 0;
 #endif
 cIndex_Set tB2 = Base2;
 for( Index h ; ( h = *(tB2++) ) < InINF ; ) {
  #if( TWOSIDED )
   if( GS[ h ] & kIsUB )
    bounds[ h ] = ub[ h ];
   else
    bounds[ h ] = lb[ h ];

   if( ABS( bounds[ h ] ) == Inf< LMNum >() )  // meanwhile, check if some
    MvdVars[ Moved++ ] = h;  // variable has in fact become unconstrained
   else
  #endif
    bNorm += bounds[ h ] * bounds[ h ];
  }

 #if( TWOSIDED )
  if( Moved ) {                       // if it has ...
   MvdVars[ Moved ] = InINF;
   CutOffConstrs( MvdVars , Moved );  // remove them
   }
 #endif

 // now compute the new Alfa[]- - - - - - - - - - - - - - - - - - - - - - - -

 RecomputeRealAlfa();

 Bf = HpINF;

 #if( LOG_BMQ > 3 )
  *BMQLog << "Modified bounds." << endl;
 #endif

 }  // end( ChangeBounds )

/*--------------------------------------------------------------------------*/

void BMinQuad::InitialSetUp( void )
{
 // set Base2[] vars- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bNorm = 0;
 const char *tGS = GS;
 Index_Set tMV = MvdVars;
 Index i = B2Dim = MB2Dim = 0;
 for( ; i < SpaceDim ; i++ ) {
  const char GSi = *(tGS++);
  if( ( GSi & kIsIn ) && ( GSi & kInB2 ) ) {
   Base2[ B2Dim++ ] = i;

   #if( TWOSIDED )
    if( GSi & kIsUB )
     bounds[ i ] = ub[ i ];
    else
     bounds[ i ] = lb[ i ];
   #endif

   cHpNum lbi = bounds[ i ];
   if( ABS( lbi ) > HpEps ) {
    *(tMV++) = i;
    bNorm += lbi * lbi;
    }
   }
  }

 Base2[ B2Dim ] = InINF;

 // now compute the Alfa[]- - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( Alfa , RealAlfa , NxtBIdx );

 if( tMV > MvdVars ) {
  *tMV = InINF;
  GiTLB( Alfa , bounds , MvdVars , Index( tMV - MvdVars ) , false );
  }

 // set MBase2[] vars - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MBase2 = Base2 + SpaceDim + 1;
 *MBase2 = InINF;

 for( ; i-- ; ) {
  const char GSi = *(--tGS);
  if( ( GSi & kIsIn ) && ( ! ( GSi & kInB2 ) ) ) {
   *(--MBase2) = i;
   MB2Dim++;
   }
  }

 // set NNStop- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( tGS = GS + ( i = SpaceDim ) ; i ; i-- )
  if( ( *(--tGS) & kENN ) == kENN )
   break;

 NNStop = i;

 // final operations- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 B2HasChgd();

 if( ActBDim ) {  // if there are items in the Bundle
  AlfaChanged();  // signal that Alfa[] has changed
  ChangeQ();      // update Q in response to the change
  UpdateB();      // update L in response to the change
  }

 Bf = HpINF;

 #if( LOG_BMQ > 3 )
  *BMQLog << "Created " << MB2Dim + B2Dim << " variables;" << endl
          << "Initial base of " << B2Dim << " variables built." << endl;
 #endif

 CheckDS();

 }  // end( InitialSetUp )

/*--------------------------------------------------------------------------*/

void BMinQuad::AddVars( cIndex_Set whch , cIndex hwmny )
{
 assert( whch[ hwmny ] == InINF );

 // count how many new variables go to [M]Base2[] - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index i;
 Index tB2Dim = 0;
 cIndex_Set tw = whch;
 for( ; ( i = *(tw++) ) < InINF ; ) {
  assert( ! ( GS[ i ] & kIsIn ) );
  GS[ i ] |= kIsIn;

  if( GS[ i ] & kInB2 )
   tB2Dim++;
  }

 cIndex tMB2Dim = hwmny - tB2Dim;

 // update NNStop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( tw = whch + hwmny ; tw > whch ; )
  if( ( GS[ *(--tw) ] & kENN ) == kENN ) {
   if( *tw >= NNStop )
    NNStop = *tw + 1;

   break;
   }

 // if necessary, construct the lists of variables in MvdVars[]- - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( tB2Dim && tMB2Dim ) {
  Index_Set tMV1 = MvdVars;
  Index_Set tMV2 = MvdVars + tB2Dim + 1;
  for( tw = whch ; ( i = *(tw++) ) < InINF ; )
   if( GS[ i ] & kInB2 )
    *(tMV1++) = i;
   else
    *(tMV2++) = i;

  *tMV1 = InINF;
  *tMV2 = InINF;
  }

 if( tB2Dim )  // if Base2 has to be updated - - - - - - - - - - - - - - - - -
 {             //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // update Base2[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set MVD = tMB2Dim ? MvdVars : whch;

  AddToB2( MVD , tB2Dim );

  // select the basic variables with nonzero bounds- - - - - - - - - - - - - -
  // meanwhile, update bounds[] and bNorm if necessary

  Index_Set tMV = MvdVars;
  for( ; ( i = *(MVD++) ) < InINF ; ) {
   #if( TWOSIDED )
    if( GS[ i ] & kIsUB )
     bounds[ i ] = ub[ i ];
    else
     bounds[ i ] = lb[ i ];
   #endif

   cHpNum bi = bounds[ i ];
   if( ABS( bi ) > HpEps ) {
    *(tMV++) = i;
    #if( ! BEXACT )
     bNorm += bi * bi;
    #endif
    }
   }

  // update Alfa[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( tMV > MvdVars ) {
   *tMV = InINF;
   GiTLB( Alfa , bounds , MvdVars , Index( tMV - MvdVars ) , false );

   AlfaChanged();
   }
  }  // end( if( tB2Dim ) )

 if( tMB2Dim )  // if MBase2 has to be updated- - - - - - - - - - - - - - - -
 {              //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // update MBase2[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set MVD = tB2Dim ? MvdVars + tB2Dim + 1 : whch;

  AddToMB2( MVD , tMB2Dim );

  // update Q and, if convenient, the base- - - - - - - - - - - - - - - - - -

  const bool UptdB = ( tMB2Dim < 2 * ReadBDim() / 3 );
  // true if updating the base is convenient

  for( ; ( i = *(MVD++) ) < InINF ; )
   AddSGSpaceDim( GiTilde( i ) , UptdB );

  if( ! UptdB )  // the base has not been updated
   UpdateB();    // do it now

  }  // end( if( tMB2Dim ) )

 B2HasChgd();
 Bf = HpINF;

 #if( LOG_BMQ > 3 )
  *BMQLog << "Created " << hwmny << " new variables" << endl;
 #endif

 CheckDS();

 }  // end( AddVars( set ) )

/*--------------------------------------------------------------------------*/

void BMinQuad::AddVars( cIndex strt , cIndex hwmny )
{
 // in this version, all the names of the variables must be larger than any
 // name of a currently defined variable

 assert( ( ! B2Dim ) || ( Base2[ B2Dim - 1 ] < strt ) );
 assert( ( ! MB2Dim ) || ( MBase2[ MB2Dim - 1 ] < strt ) );

 // first pass- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // count how many new variables go to [M]Base2[], construct the list of- - -
 // basic variables with nonzero bounds in MvdVars[], if necessary
 // update bounds[] and bNorm

 Index i;
 Index tB2Dim = 0;
 Index_Set tMV = MvdVars;
 for( i = strt ; i < strt + hwmny ; i++ ) {
  assert( ! ( GS[ i ] & kIsIn ) );
  GS[ i ] |= kIsIn;

  if( GS[ i ] & kInB2 ) {
   tB2Dim++;
   #if( TWOSIDED )
    if( GS[ i ] & kIsUB )
     bounds[ i ] = ub[ i ];
    else
     bounds[ i ] = lb[ i ];
   #endif

   cHpNum bi = bounds[ i ];
   if( ABS( bi ) > HpEps ) {
    #if( ! BEXACT )
     bNorm += bi * bi;
    #endif
    *(tMV++) = i;
    }
   }
  }

 cIndex tMB2Dim = hwmny - tB2Dim;

 // update NNStop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( i = strt + hwmny ; i-- > strt ; )
  if( ( GS[ i ] & kENN ) == kENN ) {
   if( i >= NNStop )
    NNStop = i + 1;

   break;
   }

 // update Alfa[] - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( tMV > MvdVars ) {
  *tMV = InINF;
  GiTLB( Alfa , bounds , MvdVars , Index( tMV - MvdVars ) , false );

  AlfaChanged();
  }

 // update Base2[] and MBase2[] - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // meanwhile, update Q and, if convenient, the base

 const bool UptdB = ( tMB2Dim < 2 * ReadBDim() / 3 );
 // true if updating the base is convenient

 MBase2 -= tMB2Dim;
 ShiftVect( MBase2 , MB2Dim , tMB2Dim );
 Index_Set tB2w = Base2 + B2Dim;
 Index_Set tMB2w = MBase2 + MB2Dim;
 for( i = strt ; i < strt + hwmny ; i++ )
  if( GS[ i ] & kInB2 )
   *(tB2w++) = i;
  else {
   AddSGSpaceDim( GiTilde( i ) , UptdB );
   *(tMB2w++) = i;
   }

 *tB2w = InINF;

 B2Dim += tB2Dim;
 MB2Dim += tMB2Dim;

 if( ! UptdB )  // the base has not been updated
  UpdateB();    // do it now

 B2HasChgd();
 Bf = HpINF;

 #if( LOG_BMQ > 3 )
  *BMQLog << "Created " << hwmny << " new variables from " << strt << endl;
 #endif

 CheckDS();

 }  // end( AddVars( range ) )

/*--------------------------------------------------------------------------*/

void BMinQuad::RemoveVars( cIndex_Set whch , Index hwmny )
{
 assert( whch[ hwmny ] == InINF );

 // count how many variables are removed from [M]Base2[]- - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index i;
 Index tB2Dim = 0;
 cIndex_Set tw = whch;
 for( ; ( i = *(tw++) ) < InINF ; ) {
  assert( GS[ i ] & kIsIn );
  if( GS[ i ] & kInB2 )
   tB2Dim++;
  }

 cIndex tMB2Dim = hwmny - tB2Dim;

 // if necessary, construct the lists of variables in MvdVars[]- - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( tB2Dim && tMB2Dim ) {
  Index_Set tMV1 = MvdVars;
  Index_Set tMV2 = MvdVars + tB2Dim + 1;
  for( tw = whch ; ( i = *(tw++) ) < InINF ; ) {
   if( GS[ i ] & kInB2 )
    *(tMV1++) = i;
   else
    *(tMV2++) = i;

   GS[ i ] &= kIsNN;  // keep the first bit
   }

  *tMV1 = InINF;
  *tMV2 = InINF;
  }
 else
  for( tw = whch ; ( i = *(tw++) ) < InINF ; )  // declare them all out
   GS[ i ] &= kIsNN;                            // keep the first bit

 // update NNStop- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 while( NNStop && ( ( GS[ NNStop - 1 ] & kENN ) != kENN ) )
  NNStop--;

 if( tB2Dim )  // if Base2 has to be updated - - - - - - - - - - - - - - - - -
 {             //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // update Base2[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set MVD = tMB2Dim ? MvdVars : whch;

  RmvFrmB2( MVD , tB2Dim );

  // select the basic variables with nonzero bounds- - - - - - - - - - - - - -
  // meanwhile, update bNorm if necessary

  Index_Set tMV = MvdVars;
  for( ; ( i = *(MVD++) ) < InINF ; ) {
   cHpNum bi = bounds[ i ];
   if( ABS( bi ) > HpEps ) {
    *(tMV++) = i;
    #if( ! BEXACT )
     bNorm -= bi * bi;
    #endif
    }
   }

  // update Alfa[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cIndex MVd = tMV - MvdVars;  // n. of basic variables with nonzero bound
  if( MVd ) {
   if( MVd < B2Dim ) {  // updating Alfa[] is convenient
    *tMV = InINF;
    GiTLB( Alfa , bounds , MvdVars , MVd , true );
    AlfaChanged();
    }
   else                 // recomputing Alfa[] from scratch is convenient
    RecomputeRealAlfa();

   }
  }  // end( if( tB2Dim ) )

 if( tMB2Dim )  // if MBase2 has to be updated- - - - - - - - - - - - - - - -
 {              //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // update MBase2[]- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cIndex_Set MVD = tB2Dim ? MvdVars + tB2Dim + 1 : whch;

  RmvFrmMB2( MVD , tMB2Dim );

  if( tMB2Dim < MB2Dim ) {  // if convenient, update Q- - - - - - - - - - - -
   const bool UptdB = ( tMB2Dim < 2 * ReadBDim() / 3 );
   // true if updating the base is convenient

   for( ; ( i = *(MVD++) ) < InINF ; )
    CutSGSpaceDim( GiTilde( i ) , UptdB );

   if( ! UptdB )  // the base has not been updated
    UpdateB();    // do it now
   }
  else {  // recompute Q from scratch (it costs less)- - - - - - - - - - - - -
   ChangeQ();

   if( MB2Dim + 2 < ReadBDim() )  // a clear sign that the old base is bad
    ResetB();                     // restart from scratch
   else
    UpdateB();                    // update the base

   }  // end( else( recompute Q from scratch ) ) - - - - - - - - - - - - - - -
  }  // end( if( tMB2Dim ) )

 B2HasChgd();
 Bf = HpINF;

 #if( LOG_BMQ > 3 )
  *BMQLog << "Eliminated " << hwmny << " variables" << endl;
 #endif

 CheckDS();

 }  // end( RemoveVar )

/*--------------------------------------------------------------------------*/

void BMinQuad::MoveVar( cIndex i , cIndex j , const bool iIsLst )
{
 if( ( GS[ j ] & kIsIn ) || ( ! ( GS[ i ] & kIsIn ) ) )
  return;

 if( GS[ i ] & kInB2 ) {
  Index h = iIsLst ? B2Dim - 1 : BinSearch1( Base2 , B2Dim , i );
  Index k = BinSearch2( Base2 , B2Dim , j );

  if( h < k )
   ShiftVect( Base2 + h , (--k) - h );
  else
   ShiftRVect( Base2 + k , h - k );

  Base2[ k ] = j;
  #if( TWOSIDED )
   bounds[ j ] = bounds[ i ];
  #endif
  }
 else {
  Index h = iIsLst ? MB2Dim - 1 : BinSearch1( MBase2 , MB2Dim , i );
  Index k = BinSearch2( MBase2 , MB2Dim , j );

  if( h < k )
   ShiftVect( MBase2 + h , (--k) - h );
  else
   ShiftRVect( MBase2 + k , h - k );

  MBase2[ k ] = j;
  }

 di[ j ] = di[ i ];
 lb( j ) = lb( i );
 #if( TWOSIDED )
  ub[ j ] = ub[ i ];
 #endif
 char tgsj = GS[ j ] & kIsNN;  // keep the first bit of j
 GS[ j ] = GS[ i ];
 GS[ i ] = tgsj;

 if( GS[ j ] & kIsNN ) {
  if( j >= NNStop )
   NNStop = j + 1;
  else
   while( NNStop && ( ( GS[ NNStop - 1 ] & kENN ) != kENN ) )
    NNStop--;
  }

 B2HasChgd();

 #if( LOG_BMQ > 3 )
  *BMQLog << "Variable " << i << " is now renamed as " << j << "." << endl;
 #endif

 }  // end( MoveVar )

/*--------------------------------------------------------------------------*/

void BMinQuad::RenameVars( cIndex_Set whch )
{
 if( *whch == InINF )  // no variables to remove
  return;              // nothing to do

 Index i = *(whch++);  // first variable to remove

 // find the position in Base2[] and MBase2[] - - - - - - - - - - - - - - - -

 Index_Set tB2 = Base2;
 while( *tB2 < i )
  tB2++;

 assert( *tB2 > i );

 Index_Set tMB2 = MBase2;
 while( *tMB2 < i )
  tMB2++;

 assert( *tMB2 > i );

 assert( ! ( GS[ i ] & kIsIn ) );

 // rename the variables- - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index j = i ; ++j < SpaceDim ; )
  if( *whch == j ) {
   assert( ! ( GS[ j ] & kIsIn ) );
   whch++;
   }
  else {
   const char GSi = GS[ j ];
   if( GSi & kIsIn ) {
    if( GSi & kInB2 )
     *(tB2++) = i;
    else
     *(tMB2++) = i;

    di[ i ] = di[ j ];
    #if( TWOSIDED )
     lb[ i ] = lb[ j ];
     ub[ i ] = ub[ j ];
    #endif
    bounds[ i ] = bounds[ j ];
    }

   GS[ i++ ] = GSi;
   }

 assert( *tB2 == InINF );
 assert( *tMB2 == InINF );

 // now clean the bottom variables- - - - - - - - - - - - - - - - - - - - - -

 for( ; i < SpaceDim ; i++ ) {
  GS[ i ] = kNNN;
  di[ i ] = 0;
  lb( i ) = 0;
  #if( TWOSIDED )
   ub[ i ] = Inf< LMNum >();
  #endif
  }

 // update NNStop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 while( NNStop && ( ( GS[ NNStop - 1 ] & kENN ) != kENN ) )
  NNStop--;

 B2HasChgd();

 #if( LOG_BMQ > 3 )
  *BMQLog << "Variables renamed." << endl;
 #endif

 CheckDS();

 }  // end( RenameVars )

/*--------------------------------------------------------------------------*/

void BMinQuad::MakeVarCnstr( cIndex i )
{
 GS[ i ] |= kIsNN;
 if( i >= NNStop )
  NNStop = i + 1;

 if( GS[ i ] & kInB2 ) {
  Index h = BinSearch1( MBase2 , MB2Dim-- , i );
  ShiftRVect( MBase2++ , h );

  #if( TWOSIDED )
   if( GS[ i ] & kIsUB )
    bounds[ i ] = ub[ i ];
   else
    bounds[ i ] = lb[ i ];
  #endif

  #if( ! BEXACT )
   bNorm += bounds[ i ] * bounds[ i ];
  #endif

  CutSGSpaceDim( GiTilde( i ) , bounds[ i ] );

  h = BinSearch2( Base2 , B2Dim , i );
  ShiftRVect( Base2 + h , (++B2Dim) - h );
  Base2[ h ] = i;

  B2HasChgd();
  Bf = HpINF;
  }

 #if( LOG_BMQ > 3 )
  *BMQLog << "Added constraint ";

  #if( TWOSIDED )
   if( lb[ i ] > - Inf< LMNum >() )
    *BMQLog << "[" << lb[ i ] << ", ";
   else
    *BMQLog << "[-INF, ";

   if( ub[ i ] < Inf< LMNum >() )
    *BMQLog << ub[ i ] << "]";
   else
    *BMQLog << "+INF]";
  #else
   *BMQLog << " >= " << bounds[ i ];
  #endif

  *BMQLog << " on variable " << i << endl;
 #endif

 }  // end( MakeVarCnstr )

/*--------------------------------------------------------------------------*/

void BMinQuad::MakeVarUnCnstr( cIndex i )
{
 if( GS[ i ] & kInB2 ) {
  Index h = BinSearch1( Base2 , B2Dim , i );
  ShiftVect( Base2 + h , (B2Dim--) - h );

  h = BinSearch2( MBase2 , MB2Dim++ , i );
  ShiftVect( --MBase2 , h );
  MBase2[ h ] = i;

  AddSGSpaceDim( GiTilde( i ) , bounds[ i ] );

  GS[ i ] = kEUC;

  B2HasChgd();

  #if( ! BEXACT )
   bNorm -= bounds[ i ] * bounds[ i ];
  #endif
  }
 else
  GS[ i ] &= ~kIsNN;

 while( NNStop && ( ( GS[ NNStop - 1 ] & kENN ) != kENN ) )
  NNStop--;

 #if( LOG_BMQ > 3 )
  *BMQLog << "Variable " << i << " is now unconstrained." << endl;
 #endif

 }  // end( MakeVarUnCnstr )

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SOLVING THE PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

MinQuad::MQError BMinQuad::SolveQP( HpNum ti )
{
 if( BMQt )
  BMQt->Start();

 HpNum LastBf = HpINF;
 bool AugEps = false;

 for(;;) {  // problems handling loop - - - - - - - - - - - - - - - - - - - -
  CalcOptDir( ti );

  if( ( ! QPStatus ) || ( QPStatus == kFatal ) || 
      ( QPStatus == kQPPrimUnbndd ) )
   break;

  ClearTabooList();

  if( ( Bf < HpINF ) &&
   //!! ( LastBf > Bf + eD * ( ABS( Bf ) / SpaceDim ) ) ) {
      ( LastBf > Bf + eD * std::max( double( 1 ) , std::abs( Bf ) ) / SpaceDim
	) ) {
    LastBf = Bf;     // reset the counters if there have been
    AugEps = false;  // at least one strictly decreasing step
    }

  if( ! AugEps ) {   // the first time (after a decrease)
   ChangeQ();        // refresh the data structures for Q
   UpdateB();        // therefore, refresh those for L also
   AugEps = true;    // signal that it has already been done
   continue;         // and retry;
   }

  eD *= MFactor;     // eD increase (decrease precision)
  #if( LOG_BMQ )
   *BMQLog << "Increasing eD to " << eD << endl;
  #endif

  if( eD > MaxEpsD ) { 
   #if( LOG_BMQ )
    *BMQLog << endl << "ERROR: eD (" << eD << ") out of limits. [Solve[B]QP]" 
            << endl;
   #endif

   QPStatus = kFatal;
   break;
   }
  }  // end while( problems ) - - - - - - - - - - - - - - - - - - - - - - - -

 if( BMQt )
  BMQt->Stop();

 return( QPStatus );
 
 }  // end( BMinQuad::SolveQP )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR READING RESULTS ---------------------*/
/*--------------------------------------------------------------------------*/

void BMinQuad::ReadZ( LMRow g )
{
 cIndex_Set o = Base2;
 Index h;

 for( ; ( h = *(o++) ) < InINF ; )
  g[ h ] = di[ h ];

 for( o = MBase2 ; ( h = *(o++) ) < InINF ; )
  g[ h ] = di[ h ];

 }  // end( ReadZ )

/*--------------------------------------------------------------------------*/

void BMinQuad::ReadZSprs( LMRow g )
{
 cIndex_Set tB2 = Base2;
 cIndex_Set tMB2 = MBase2;
 Index h = *tB2;
 Index k = *tMB2;

 while( h != k )
  if( h < k ) {
   *(g++) = di[ h ];
   h = *(++tB2);
   }
  else {
   *(g++) = di[ k ];
   k = *(++tMB2);
   }
 }  // end( ReadZSprs )

/*--------------------------------------------------------------------------*/

Index BMinQuad::ReadVNames( Index_Set VNames )
{
 cIndex_Set tB2 = Base2;
 cIndex_Set tMB2 = MBase2;
 Index h = *tB2;
 Index k = *tMB2;

 while( h != k )
  if( h < k ) {
   *(VNames++) = h;
   h = *(++tB2);
   }
  else {
   *(VNames++) = k;
   k = *(++tMB2);
   }

 return( B2Dim + MB2Dim );

 }  // end( ReadVNames )

/*--------------------------------------------------------------------------*/

void BMinQuad::ReadD( LMRow d , cIndex CpyFrst )
{
 cIndex_Set o = Base2;
 Index h;

 if( CpyFrst ) {
  h = *(o++);
  for( Index i = 0 ; i < CpyFrst ; i++ )
   if( i == h ) {
    d[ i ] = bounds[ h ];
    h = *(o++);
    }
   else
    d[ i ] = - PrvsTi * di[ i ];
  }
 else {
  for( ; ( h = *(o++) ) < InINF ; )
   d[ h ] = bounds[ h ];

  for( o = MBase2 ; ( h = *(o++) ) < InINF ; )
   d[ h ] = - PrvsTi * di[ h ];
  }
 }  // end( ReadD )

/*--------------------------------------------------------------------------*/

void BMinQuad::ReadDSprs( LMRow d )
{
 cIndex_Set tB2 = Base2;
 cIndex_Set tMB2 = MBase2;
 Index h = *tB2;
 Index k = *tMB2;

 while( h != k )
  if( h < k ) {
   *(d++) = bounds[ h ];
   h = *(++tB2);
   }
  else {
   *(d++) = - PrvsTi * di[ k ];
   k = *(++tMB2);
   }
 }  // end( ReadDSprs )

/*--------------------------------------------------------------------------*/

HpNum BMinQuad::DPerG( cSgRow g )
{
 HpNum scpr = 0;
 cIndex_Set o = Base2;
 Index h;

 for( ; ( h = *(o++) ) < InINF ; )
  scpr -= g[ h ] * bounds[ h ];

 scpr /= PrvsTi;

 for( o = MBase2 ; ( h = *(o++) ) < InINF ; )
  scpr += g[ h ] * di[ h ];

 return( scpr );

 }  // end( DPerG )

/*--------------------------------------------------------------------------*/

void BMinQuad::AddD( LMRow L1 , cLMRow L2 , cHpNum Tau )
{
 cIndex_Set o = Base2;
 Index h;

 if( Tau == PrvsTi )
  for( ; ( h = *(o++) ) < InINF ; )
   L1[ h ] = L2[ h ] + bounds[ h ];
 else
  for( HpNum TauR = Tau / PrvsTi ; ( h = *(o++) ) < InINF ; )
   L1[ h ] = L2[ h ] + TauR * bounds[ h ];

 for( o = MBase2 ; ( h = *(o++) ) < InINF ; )
  L1[ h ] = L2[ h ] - Tau * di[ h ];

 }  // end( AddD )

/*--------------------------------------------------------------------------*/

void BMinQuad::AddDSprs( LMRow L1 , cLMRow L2 , cHpNum Tau )
{
 cIndex_Set tB2 = Base2;
 cIndex_Set tMB2 = MBase2;
 Index h = *tB2;
 Index k = *tMB2;

 if( Tau == PrvsTi ) {
  while( h != k )
   if( h < k ) {
    *(L1++) = *(L2++) + bounds[ h ];
    h = *(++tB2);
    }
   else {
    *(L1++) = *(L2++) - Tau * di[ k ];
    k = *(++tMB2);
    }
   }
 else {
  HpNum TauR = Tau / PrvsTi;

  while( h != k )
   if( h < k ) {
    *(L1++) = *(L2++) + TauR * bounds[ h ];
    h = *(++tB2);
    }
   else {
    *(L1++) = *(L2++) - Tau * di[ k ];
    k = *(++tMB2);
    }
  }
 }  // end( AddDSprs )

/*--------------------------------------------------------------------------*/

HpNum BMinQuad::ReadGTz( cIndex i )
{
 return( MinQuad::ReadGTz( i ) - ( RealAlfa[ i ] - Alfa[ i ] ) / PrvsTi );
 }

/*--------------------------------------------------------------------------*/

void BMinQuad::SensitAnals1( HpNum &v1 , HpNum &v2 , HpNum &v3 )
{
 MinQuad::SensitAnals1( v1 , v2 , v3 );
 v3 += bNorm;
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

BMinQuad::~BMinQuad( void )
{
 // First: print out information - - - - - - - - - - - - - - - - - - - - - -

 #if( LOG_BMQ > 2 )
  *BMQLog << endl << "BCalls " << BCalls << " (faults " << BCalls - BSccss
          << ") ~ " << SpaceDim << " variables, " << SumAverages / BCalls
	  << " avg. constr.";
  if( BMQt )
   *BMQLog << ", " << BMQt->Read() << "time";

  *BMQLog << endl;
 #elif( LOG_BMQ > 1 )
  *BMQLog << SpaceDim << "\t" << BCalls << "\t" << SumAverages / BCalls
	  << "\t";
  if( BMQt )
   *BMQLog  << BMQt->Read() << "\t";
 #endif

 // Second: memory deallocation - - - - - - - - - - - - - - - - - - - - - - -

 delete BMQt;

 if( MaxBDim )
  MemDealloc();

 }  // end( ~BMinQuad )

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

void BMinQuad::RecomputeRealAlfa( void )
{
 VectAssign( Alfa , RealAlfa , NxtBIdx );
 GiTLB( Alfa , bounds , Base2 , B2Dim , false );

 AlfaChanged();
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

inline void BMinQuad::CalcOptDir( HpNum ti )
{
 /* Same as in the base class, but taking into account the box constraints
    (in a "smart" way): calls MinQuad::SolveQP() as a subroutine. */

 #if( LOG_BMQ > 1 )
  BCalls++;
  float TmpSumAverages = 0;
  unsigned long int BStep = 0;
 #endif

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /*- - - - - - - - - - - - - - Main Loop - - - - - - - - - - - - - - - - - -*/
 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 HpNum BfBot = HpINF;  // "security level", used to switch between the "quick
                       // and dirty" pricing and the "slow but clean" one
 HpNum MnDcrs = 0;     // minimum decrease to declare a descent step

 Index tmpMVA = std::min( MaxVarAdd , B2Dim + MB2Dim );
                       // max number of variables that can be priced in
                       // in the same iteration; in case of problems it is
                       // gradually decreased down to one ("strict pricing")
                       // to avoid loops (together with the "taboo" mechanism)
 Index Moved = 0;      // number of variables moved from Base2 to MBase2 or
                       // vice-versa
 Index Entrd = InINF;  // name of a variable entered in Base2 in the
                       // latest iteration
 do {
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - - Inner Loop  - - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Within this inner loop, non-violated box constraints are removed (in
    other words, negative Xsi[]'s are eliminated).
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  do {
   #if( LOG_BMQ > 1 )
    #if( LOG_BMQ > 3 )
     *BMQLog << endl << "[" << BCalls << "/" << BStep++ << "-" << B2Dim
	     << "+" << MB2Dim << "]: bNorm = " << bNorm << " ~ Bf = ";
     if( Bf < HpINF ) *BMQLog << Bf << endl;
     else             *BMQLog << "INF" << endl;
    #endif

    TmpSumAverages += B2Dim;
   #endif

   CheckDS();

   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // the "usual" (QP) is solved- - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   MinQuad::SolveQP( ti );

   if( ( QPStatus == kFatal ) ||
       ( ( QPStatus == kQPPrimUnbndd ) && ( ! B2Dim ) ) ) {
    #if( LOG_BMQ > 1 )
     SumAverages += TmpSumAverages / BStep;
    #endif
    return;  // errors are just passed up
    }

   if( QPStatus == kQPPrimUnbndd ) {
    // here a rather funny situation occurs: the problem has no feasible
    // primal solution (i.e., it is dual unbounded) because the variables
    // fixed to bounds in Base2[] kill every possible solution; however, a
    // different fixing (with less variables) may allow a solution to exist
    // of course this is only possible if Base2[] is nonempty

    #if( LOG_BMQ > 3 )
     *BMQLog << " Inner QP unfeasible";
    #endif

    if( Entrd == InINF ) {
     // the base was feasible at the previous call, but it has become
     // unfeasible due to "external" intervention (e.g. constraint added):
     // reset Base2[] to the empty set

     #if( LOG_BMQ > 3 )
      *BMQLog << ", resetting Base2[]" << endl;
     #endif

     Moved = B2Dim;
     B2Dim = MB2Dim = 0;
     *Base2 = InINF;
     MBase2 = Base2 + SpaceDim + 1;
     *MBase2 = InINF;

     for( char *tGS = GS + SpaceDim ; tGS-- > GS ; )
      if( *tGS & kIsIn ) {
       #if( TWOSIDED )
        *tGS &= ~ ( kInB2 | kIsUB );
       #else
        *tGS &= ~ kInB2;
       #endif
       *(--MBase2) = tGS - GS;
       MB2Dim++;
       }

     bNorm = 0;
     VectAssign( Alfa , RealAlfa , NxtBIdx );

     B2HasChgd();

     if( ActBDim ) {  // if there are items in the Bundle
      AlfaChanged();  // signal that Alfa[] has changed
      ChangeQ();      // update Q in response to the change
      UpdateB();      // update L in response to the change
      }

     continue;
     }

    // else, for the time being ignore the fact that the inner QP is
    // unfeasible and pretend that Mult[] represent a feasible d[],
    // which it doesn't; however, if the so-computed d[] is already
    // bound unfeasible, then Base2[] has to be reduced anyway
    /*!!
    else {
     // this is the result of putting some new variable in Base2: undo
     // last move and be more conservative for the future

     #if( LOG_BMQ > 3 )
      *BMQLog << ", undoing last move" << endl;
     #endif

     // MvdVars[] is already filled with the variables just put in

     #if( ! BEXACT )
      bNorm -= Norm( bounds , MvdVars );
     #endif

     VectAssign( GS , kENN , MvdVars );
     CutOffConstrs( MvdVars , Moved );

     if( tmpMVA > 1 ) {  // multiple price-in was used
      tmpMVA /= 2;       // restrict the pricing window (downto 1)
      #if( LOG_BMQ > 3 )
       *BMQLog << " setting MVA = " << tmpMVA << endl;
      #endif
      }
     else {              // this happened under single-variable price-in
      GS[ Entrd ] |= kIsTb;
      TLDim++;

      #if( LOG_BMQ > 3 )
       *BMQLog << " marking " << Entrd << " as taboo" << endl;
      #endif
      }
     continue;
     }
     !!*/
    }

   MBHasChgd();

   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // feasibility check- - - - - -- - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Moved = 0;

   if( BfBot < HpINF )  // "secure" pricing - - - - - - - - - - - - - - - - -
   {                    //- - - - - - - - - - - - - - - - - - - - - - - - - -
    #if( LAZY_D == 0 )
     CalculateZ( tmpdi );
    #elif( LAZY_D == 1 )
     CalculateZ( Base2 , tmpdi );
    #endif

    // check if tmpdi is feasible: if not, find the - - - - - - - - - - - - -
    // maximum feasible step along tmpd - di

    LMNum step = Inf< LMNum >();
    cHpNum eDir = eD * std::max( BDim , Index( 1 ) );

    cIndex_Set tB2 = Base2;
    for( Index h ; ( h = *(tB2++) ) < InINF ; ) {
     #if( LAZY_D == 2 )
      tmpdi[ h ] = CalculateZ( h );
     #endif

     cLMNum d1h = - ti * tmpdi[ h ];

     #if( TWOSIDED )
      if( GS[ h ] & kIsUB ) {
       if( d1h + size( d1h ) > ub[ h ] )
	continue;
       }
      else
     #endif
       if( d1h - size( d1h ) < lb( h ) )
	continue;

     cLMNum dh = ti * di[ h ];
     step = std::min( step , ( dh + bounds[ h ] ) / ( dh + d1h ) );
     }

    if( step == Inf< LMNum >() )  // tmpdi is feasible- - - - - - - - - - - -
     Swap( di , tmpdi );
    else {               // move along tmpd - di of step- - - - - - - - - - -
     #if( LOG_BMQ > 3 )
      *BMQLog << " max step = " << step << endl;
     #endif

     cIndex_Set tB2 = Base2;
     for( Index h ; ( h = *(tB2++) ) < InINF ; ) {
      di[ h ] += step * ( tmpdi[ h ] - di[ h ] );

      if( ABS( ti * di[ h ] + bounds[ h ] ) <= size( bounds[ h ] ) ) {
       #if( LOG_BMQ > 4 )
        *BMQLog << " Out constraint, d[ " << h << " ] = " << - ti * di[ h ]
		<< endl;
       #endif
       #if( ! BEXACT )
        bNorm -= bounds[ h ] * bounds[ h ];
       #endif

       GS[ h ] = kENN;
       MvdVars[ Moved++ ] = h;
       }
      }  // end for( tB2 )

     if( ! Moved ) {  // haven't found any constraint to eliminate- - - - - -
      #if( LOG_BMQ )
       *BMQLog << endl << "Fault: unable to reach the boundary. [BCalcOptDir]"
	       << endl;
      #endif

      QPStatus = kLoop;
      return;
      }
     }   // end else( move along tmpdi - di )
    }    // end if( "secure" pricing )
   else  // "quick and dirty" pricing - - - - - - - - - - - - - - - - - - - -
   {     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #if( LAZY_D == 0 )
     CalculateZ( di );
    #elif( LAZY_D == 1 )
     CalculateZ( Base2 , di );
    #endif

    cHpNum eDir = eD * std::max( BDim , Index( 1 ) );
    cIndex_Set tB2 = Base2;

    for( Index h ; ( h = *(tB2++) ) < InINF ; ) {
     #if( LAZY_D == 2 )
      di[ h ] = CalculateZ( h );
     #endif

     cLMNum dh = - ti * di[ h ];

     #if( TWOSIDED )
      if( GS[ h ] & kIsUB ) {
       if( dh + size( dh ) > ub[ h ] )
        continue;
       }
      else
     #endif
       if( dh - size( dh ) < lb( h ) )
        continue;

     #if( LOG_BMQ > 4 )
      #if( TWOSIDED )
       if( GS[ h ] & kIsUB )
        *BMQLog << " Out UB constraint, d[ " << h << " ] = "
		<< dh << " < " << ub[ h ] << endl;
       else
      #endif
        *BMQLog << " Out LB constraint, d[ " << h << " ] = "
		<< dh << " > " << lb( h ) << endl;
     #endif
     #if( ! BEXACT )
      bNorm -= bounds[ h ] * bounds[ h ];
     #endif

     GS[ h ] = kENN;
     MvdVars[ Moved++ ] = h;
     if( Moved >= MaxVarRmv )
      break;
     }
    }  // end else( "quick and dirty" pricing ) - - - - - - - - - - - - - - -
       // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   if( ( QPStatus == kQPPrimUnbndd ) && ( ! Moved ) ) {
    // dir[] is feasible, but the inner QP is unfeasible; that is only one
    // of infinitely many possible directions
    // check if dir[] + \beta * tmpdi[] is feasible for every possible value
    // of \beta > 0; if not, find minimum \beta which brings it on the
    // boundary (thus identifying at least one constraint to be removed)

    Swap( Mult , tmpv );  // dirty trick: this way CalculateZ() computes the
                          // variable part of dir[]
    #if( LAZY_D == 0 )
     CalculateZ( tmpdi );
    #elif( LAZY_D == 1 )
     CalculateZ( Base2 , tmpdi );
    #endif

    LMNum beta = Inf< LMNum >();
    cHpNum eDir = eD * std::max( BDim , Index( 1 ) );

    cIndex_Set tB2 = Base2;
    for( Index h ; ( h = *(tB2++) ) < InINF ; ) {
     #if( LAZY_D == 2 )
      tmpdi[ h ] = CalculateZ( h );
     #endif

     #if( TWOSIDED )
      if( GS[ h ] & kIsUB ) {
       // - t * dir[ h ] >= ub[ h ]
       // - t * ( dir[ h ] + \beta * tmpdi[ h ] ) >= ub[ h ]
       // either tmpdi[ h ] <= 0 
       // \beta <= ( ub[ h ] + t * dir[ h ] ) / ( - t * tmpdi[ h ] )
       if( tmpdi[ h ] <= eDir )
	continue;
       }
      else
     #endif
       // - t * dir[ h ] <= lb[ h ]
       // - t * ( dir[ h ] + \beta * tmpdi[ h ] ) <= lb[ h ]
       // either tmpdi[ h ] >= 0
       // \beta <= ( lb[ h ] + t * dir[ h ] ) /  ( - t * tmpdi[ h ] )
       if( tmpdi[ h ] >= - eDir )
	continue;

     LMNum betah = ( bounds[ h ] + ti * di[ h ] ) / ( - ti * tmpdi[ h ] );
     if( betah < beta )
      beta = betah;

     }  // end( for( h ) )

    Swap( Mult , tmpv );  // restore true solution

    if( beta == Inf< LMNum >() )  // the inner QP is *really* unbounded - - -
     return;
    else {               // dir = dir + \beta * tmpd- - - - - - - - - - - - -
     #if( LOG_BMQ > 3 )
      *BMQLog << " max beta = " << beta << endl;
     #endif

     cIndex_Set tB2 = Base2;
     for( Index h ; ( h = *(tB2++) ) < InINF ; ) {
      di[ h ] += beta * tmpdi[ h ];

      if( ABS( ti * di[ h ] + bounds[ h ] ) <= size( bounds[ h ] ) ) {
       #if( LOG_BMQ > 4 )
        *BMQLog << " Out constraint, d[ " << h << " ] = " << - ti * di[ h ]
		<< endl;
       #endif
       #if( ! BEXACT )
        bNorm -= bounds[ h ] * bounds[ h ];
       #endif

       GS[ h ] = kENN;
       MvdVars[ Moved++ ] = h;
       }
      }  // end for( tB2 )
 
     if( ! Moved ) {  // haven't found any constraint to eliminate
      #if( LOG_BMQ )
       *BMQLog << endl << "Fault: unable to reach the boundary. [BCalcOptDir]"
	       << endl;
      #endif

      QPStatus = kLoop;
      return;
      }
     }   // end( else( dir = dir + \beta * tmpd ) ) - - - - - - - - - - - - - 
    }  // end( if( QPStatus == kQPPrimUnbndd ... ) )

   // actually eliminate identified constraints (if any) from Base2[] - - - - 
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

   if( Moved ) {
    #if( LOG_BMQ == 4 )
     *BMQLog << " Out " << Moved << " constraints" << endl;
    #endif

    MvdVars[ Moved ] = InINF;
    CutOffConstrs( MvdVars , Moved );
    }
   #if( LOG_BMQ > 3 )
    else
     *BMQLog << " Base2 feasible";
   #endif

   } while( Moved );

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - End Inner Loop  - - - - - - - - - - - - - - -*/
  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  #if( BEXACT )
   bNorm = Norm( bounds , Base2 );
  #else
   if( ! B2Dim )
    bNorm = 0;   // from time to time, clean up bNorm
  #endif

  HpNum newBf = ( ti * Quad - bNorm / ti ) / 2 + Lin;

  #if( LOG_BMQ > 3 )
   *BMQLog << " ~ bNorm = " << bNorm << " ~ newBf = " << newBf << endl;
  #endif

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    The current base is feasible, so calculate Bf and test for decrease.
    If the "quick and dirty" rule is used, the decrease is *not guaranteed*
    and this control may fail: in this case, the "secure" rule will be used
    afterwards until a strict decrease (w.r.t. the "old" Bf at *this*
    iteration) is obtained.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  if( Bf < HpINF ) {  // skipping the tests if Bf is undefined
   // Test the relative decrease of B by checking ( Bf - newBf ) / Bf against
   // eD / SpaceDim; the idea is that each of the SpaceDim variables can make
   // Bf to decrease of at most 1 / SpaceDim of its current value

   //!! cHpNum RqrdDcrs = eD * ( ABS( Bf ) / SpaceDim );
   cHpNum RqrdDcrs = std::max( MnDcrs ,
   //!!                  HpEps * std::max( BDim , Index( 1 ) ) * ABS( Bf ) );
			       eD * std::max( double( 1 ) , std::abs( Bf ) )
			       / SpaceDim );

   if( Bf - newBf < - RqrdDcrs ) {  // an increasing step
    #if( LOG_BMQ > 3 )
     *BMQLog << " Increasing step (Df = " << newBf - Bf << "): ";
    #endif

    if( BfBot == HpINF ) {          // it's the first increasing step
     BfBot = Bf;                    // just start using the "secure" pricing

     #if( LOG_BMQ > 3 )
      *BMQLog << "starting secure pricing" << endl;
     #endif
     }
    else {                          // secure pricing is already active
     #if( LOG_BMQ > 2 )
      *BMQLog << " Fault [BCalcOptDir]" << endl;
     #endif

     QPStatus = kLoop;
     return;
     }
    }
   else
    if( Bf - newBf <= RqrdDcrs ) {  // a "short" step
     #if( LOG_BMQ > 3 )
      *BMQLog << endl << " Short step (Df = " << newBf - Bf << "),";
     #endif

     if( tmpMVA > 1 ) {    // multiple price in was used
      tmpMVA /= 2;         // restrict the pricing window (downto 1)
      #if( LOG_BMQ > 3 )
       *BMQLog << " setting MVA = " << tmpMVA << endl;
      #endif
      }
     else                  // strict pricing was already used
      if( Entrd < InINF ) {
       GS[ Entrd ] |= kIsTb;
       TLDim++;

       #if( LOG_BMQ > 3 )
        *BMQLog << " marking " << Entrd << " as taboo" << endl;
       #endif
       }
     }
    else                   // a "good step"
     //!! if( newBf + eD * ( ABS( newBf ) / SpaceDim ) < BfBot ) {
     if( newBf + eD * std::max( double( 1 ) , std::abs( newBf ) ) / SpaceDim
	 < BfBot ) {
      #if( LOG_BMQ > 3 )
       if( BfBot < HpINF )
        *BMQLog << " disabling secure pricing" << endl;
      #endif
      BfBot = HpINF;  // it is safe to disable the "secure" rule ...

      #if( LOG_BMQ > 3 )
       if( tmpMVA < std::min( MaxVarAdd , B2Dim + MB2Dim ) )
        *BMQLog << " disabling strict pricing" << endl;
      #endif
      tmpMVA = std::min( MaxVarAdd , B2Dim + MB2Dim );
      //... the "strict" pricing

      ClearTabooList();    // ... and to clear the taboo list
      }

   }  // end if( Bf < HpINF )

  Bf = newBf;  // update the function value

  /* Note: there used to be a check here to the tune that if Bf is "small
     enough", one could can detect optimality without the (costly) check of
     constraints violation. However, the entries of d[] in MBase2 would have
     to be "artificially" set to 0, which may lead to numerical errors
     afterwards. Hence, the check is better disabled.

  if( ( Bf = newBf ) <= MinFVal ) {
   #if( LOG_BMQ > 3 )
    *BMQLog << " < MinFVal = " << MinFVal << ": STOP" << endl;
   #endif

   cIndex_Set o = MBase2;
   for( Index h ; ( h = *(o++) ) < InINF ; )
    di[ h ] = 0;

   break;
   }
  */

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Now check for the existence of violated (primal) constraints to be put in
    the "constraints base" to obtain a decrease of Bf; in other words, check
    for any Xsi[ h ] == 0 that, if let free to become > 0, can diminish the
    value of Bf.
    It can be shown that putting in the "constraints base" a (violated)
    constraint h such that
      vh = lb[ h ] - di[ h ] ( di[ h ] - ub[ h ] ) > 0
    ensures a decrease in the value of Bf of at least vh^2 / ( 2 * t ).
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  MnDcrs = 0;  // lower bound on the decrease of Bf

  #if( LAZY_D == 1 )
   CalculateZ( MBase2 , di );
  #endif

  cHpNum eDir = eD * std::max( BDim , Index( 1 ) );
  cIndex_Set tMB2 = MBase2;

  for( Index h ; ( h = *tMB2 ) < NNStop ; tMB2++ )
   if( GS[ h ] == kENN ) {  // only existing NN variables in MBase2 which
    #if( LAZY_D == 2 )      // are not "taboo" can enter Base2
     di[ h ] = CalculateZ( h );
    #endif

    cLMNum dh = - ti * di[ h ];
    HpNum vh;

    if( dh + size( dh ) < lb( h ) ) {
     vh = lb( h ) - dh;
     #if( LOG_BMQ > 4 )
      *BMQLog << " In LB constraint " << h << ", vh = " << vh << endl;
     #endif

     GS[ h ] = kLBC;
     #if( TWOSIDED )
      bounds[ h ] = lb[ h ];
      }
     else
      if( dh - size( dh ) > ub[ h ] ) {
       vh = dh - ub[ h ];
       #if( LOG_BMQ > 4 )
        *BMQLog << " In UB constraint " << h << ", vh = " << vh << endl;
       #endif

       GS[ h ] = kUBC;
       bounds[ h ] = ub[ h ];
     #endif
       }
      else
       continue;

    #if( ! BEXACT )
     bNorm += bounds[ h ] * bounds[ h ];
    #endif
    if( vh > MnDcrs )
     MnDcrs = vh;

    MvdVars[ Moved++ ] = h;
    if( Moved >= tmpMVA )
     break;

    }  // end for( h )

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

  if( Moved ) {
   #if( LOG_BMQ == 4 )
    *BMQLog << " In " << Moved << " constraints" << endl;
   #endif

   Entrd = *MvdVars;       // remember the name of one of the entered items
   MnDcrs *= MnDcrs;       // compute a lower bound on the expected
   MnDcrs /= ( ti + ti );  // decrease of Bf

   MvdVars[ Moved ] = InINF;
   PutInConstrs( MvdVars , Moved );
   }

  //!!  if( ( ! Moved ) && ( Bf - eD * ( ABS( Bf ) / SpaceDim ) > newBf ) ) {
  if( ( ! Moved ) && ( Bf - eD * std::max( double( 1 ) , std::abs( Bf ) )
		       / SpaceDim > newBf ) ) {
   // if no variables enter in Base2 (some because they are "taboo"),
   // but Bf() is not (about) as good as BfBot, then kLoop is returned

   #if( LOG_BMQ > 2 )
    *BMQLog << endl << "Fault: Bf() can't reach the min. value " << BfBot
	    << " [BCalcOptDir]";
   #endif

   QPStatus = kLoop;
   return;
   }
  } while( Moved );

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /*- - - - - - - - - - - - - End Main Loop - - - - - - - - - - - - - - - - -*/
 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ClearTabooList();  // clear the taboo list

 #if( LOG_BMQ > 1 )
  SumAverages += TmpSumAverages / BStep;
  BSccss++;

  #if( LOG_BMQ > 3 )
   *BMQLog << " Stop: BOptimum reached ~ Bf = " << Bf << endl << endl;
  #endif
 #endif

 }  // end( BMinQuad::CalcOptDir )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::AddToB2( cIndex_Set MVD , cIndex MVDd )
{
 if( B2Dim ) {
  cIndex_Set tB2r = Base2 + B2Dim;
  Index_Set tB2w = Base2 + B2Dim + MVDd;
  *(tB2w--) = InINF;
  cIndex_Set tMVD = MVD + MVDd;
  Index h = *(--tMVD);

  if( *MVD < *Base2 ) {  // the new first element in Base2[] comes from MVD[]
   for( Index k = *(--tB2r) ; ; )
    if( h > k ) {
     *(tB2w--) = h;
     h = *(--tMVD);
     }
    else {
     *(tB2w--) = k;
     if( tB2r > Base2 )
      k = *(--tB2r);
     else
      break;
     }

   for( *(tB2w--) = h ; tMVD > MVD ; )
    *(tB2w--) = *(--tMVD);
   }
  else                 // the first element of Base2[] remains the same
   for( Index k = *(--tB2r) ; ; )
    if( h > k ) {
     *(tB2w--) = h;
     if( tMVD > MVD )
      h = *(--tMVD);
     else
      break;
     }
    else {
     *(tB2w--) = k;
     k = *(--tB2r);
     }
  }
 else
  VectAssign( Base2 , MVD , MVDd + 1 );

 B2Dim += MVDd;

 CheckB2();

 }  // end( AddToB2 )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::RmvFrmB2( cIndex_Set MVD , cIndex MVDd )
{
 Index_Set tB2w = Base2 + BinSearch1( Base2 , B2Dim , *MVD );
 cIndex_Set tB2r = tB2w + 1;
 Index h  = *(++MVD);

 for( ; h < InINF ; tB2r++ )
  if( *tB2r == h )
   h = *(++MVD);
  else
   *(tB2w++) = *tB2r;

 while( ( h = *(tB2r++) ) < InINF )
  *(tB2w++) = h;

 *tB2w = InINF;

 B2Dim -= MVDd;

 CheckB2();

 }  // end( RmvFrmB2 )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::AddToMB2( cIndex_Set MVD , cIndex MVDd )
{
 Index_Set tMB2w = MBase2 - MVDd;
 Index_Set tMB2r = MBase2;
 Index k = *(tMB2r++);

 for( Index h = *(MVD++) ; ; tMB2w++ )
  if( h < k ) {
   *tMB2w = h;
   if( ( h = *(MVD++) ) == InINF )
    break;
   }
  else {
   *tMB2w = k;
   k = *(tMB2r++);
   }

 MBase2 -= MVDd;
 MB2Dim += MVDd;

 CheckMB2();

 }  // end( AddToMB2 )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::RmvFrmMB2( cIndex_Set MVD , cIndex MVDd )
{
 cIndex_Set tMVD = MVD + MVDd;
 Index_Set tMB2w = MBase2 + BinSearch1( MBase2 , MB2Dim , *(--tMVD) );
 cIndex_Set tMB2r = tMB2w - 1;

 for( ; tMVD > MVD ; )
  for( cIndex h = *(--tMVD) ;; ) {
   cIndex k = *(tMB2r--);
   if( k > h )
    *(tMB2w--) = k;
   else
    break;
   }

 while( tMB2r >= MBase2 )
  *(tMB2w--) = *(tMB2r--);

 MBase2 += MVDd;
 MB2Dim -= MVDd;

 CheckMB2();

 }  // end( RmvFrmMB2 )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::CutOffConstrs( cIndex_Set MVD , Index MVDd )
{
 // removes the basic constraints on the MVDd variables in MVD[] from Base2,
 // putting them in MBase2; does not handle bNorm
 // note that it assumes that GS[] has already been updated, so at the
 // beginning of the call GS[] and Base2[]/MBase2[] are out of sync

 assert( MVDd > 0 );
 assert( MVD[ MVDd ] == InINF );

 // first eliminate them all from Base2[] - - - - - - - - - - - - - - - - - -

 RmvFrmB2( MVD , MVDd );

 // then insert them all in MBase2[]- - - - - - - - - - - - - - - - - - - - -

 AddToMB2( MVD , MVDd );

 CheckGS();

 // now take care of the rest - - - - - - - - - - - - - - - - - - - - - - - -

 const bool UptdB = ( MVDd < 2 * ReadBDim() / 3 );
 // true if updating the base is convenient

 for( Index h ; ( h = *(MVD++) ) < InINF ; )
  AddSGSpaceDim( GiTilde( h ) , bounds[ h ] , UptdB );

 if( ! UptdB )  // the base has not been updated
  UpdateB();    // do it now

 B2HasChgd();

 }  // end( CutOffConstrs )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::PutInConstrs( cIndex_Set MVD , Index MVDd )
{
 // put in Base2[] the constraints on the MVDd variables in MVD[], removing
 // them from MBase2[]; does not handle bNorm
 // note that it assumes that GS[] has already been updated, so at the
 // beginning of the call GS[] and Base2[]/MBase2[] are out of sync

 assert( MVDd > 0 );
 assert( MVD[ MVDd ] == InINF );

 // first eliminate them all from MBase2[]- - - - - - - - - - - - - - - - - -

 RmvFrmMB2( MVD , MVDd );

 // then insert them all in Base2[] - - - - - - - - - - - - - - - - - - - - -

 AddToB2( MVD , MVDd );

 CheckGS();

 // now take care of the rest- - - - - - - - - - - - - - - - - - - - - - - - -

 if( MVDd > MB2Dim ) {  // updating Q would cost more than recomputing it
  ChangeQ();            // ... so recompute it

  // now put in MvdVars[] the subset of MVD[] corresponding to variables
  // with nonzero bound; note that MVD[] may well *be* MvdVars[]

  Index_Set tMV = MvdVars;
  for( Index h ; ( h = *(MVD++) ) < InINF ; )
   if( ABS( bounds[ h ] ) > HpEps )
    *(tMV++) = h;

  if( tMV > MvdVars ) {  // if there is any
   *tMV = InINF;
   GiTLB( Alfa , bounds , MvdVars , Index( tMV - MvdVars ) , false );
   AlfaChanged();
   }

  if( MB2Dim + 2 < ReadBDim() )  // a clear sign that the old base is bad
   ResetB();                     // restart from scratch
  else
   UpdateB();                    // update the base
  }
 else {                 // updating Q is convenient
  const bool UptdB = ( MVDd < 2 * ReadBDim() / 3 );
  // true if updating the base is convenient
  for( Index h ; ( h = *(MVD++) ) < InINF ; ) 
   CutSGSpaceDim( GiTilde( h ) , bounds[ h ] , UptdB );

  if( ! UptdB )  // the base has not been updated
   UpdateB();    // do it now
  }

 CheckRA();

 B2HasChgd();

 }  // end( PutInConstrs )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::ClearTabooList( void )
{
 if( TLDim ) {
  cIndex_Set tB = MBase2;
  Index h;

  for( ; ( h = *(tB++) ) < InINF ; )
   if( GS[ h ] & kIsTb ) {
    GS[ h ] &= ~kIsTb;
    if( ! --TLDim )
     break;
    }

  if( TLDim )
   for( tB = Base2 ; ( h = *(tB++) ) < InINF ; )
    if( GS[ h ] & kIsTb ) {
     GS[ h ] &= ~kIsTb;
     if( ! --TLDim )
      break;
     }

  #if( LOG_BMQ > 3 )
   *BMQLog << " taboo list cleared" << endl;
  #endif
  }
 }  // end( ClearTabooList )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::MemDealloc( void )
{
 #if( TWOSIDED )
  delete[] ub;
  delete[] lb;
 #endif

 delete[] bounds;
 #if( ! CNDVD_TMP )
  delete[] tmpdi;
 #endif
 delete[] di;

 delete[] MvdVars;
 delete[] Base2;
 delete[] GS;

 delete[] RealAlfa;

 }  // end( MemDealloc )

/*--------------------------------------------------------------------------*/

inline void BMinQuad::CheckB2( void )
{
 #if( CHECK_DS & 1 )
  assert( Base2[ B2Dim ] == InINF );

  if( B2Dim ) {
   cIndex_Set tB2 = Base2;
   Index i = *(tB2++);
   for( Index h ; ( h = *(tB2++) ) < InINF ; ) {
    assert( i < h );
    i = h;
    }

   assert( Index( tB2 - Base2 ) == B2Dim + 1 );
   }
 #endif

 }  // end( CheckB2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline void BMinQuad::CheckMB2( void )
{
 #if( CHECK_DS & 1 )
  assert( MBase2[ MB2Dim ] == InINF );

  if( MB2Dim ) {
   cIndex_Set tMB2 = MBase2;
   Index i = *(tMB2++);
   for( Index h ; ( h = *(tMB2++) ) < InINF ; ) {
    assert( i < h );
    i = h;
    }

   assert( Index( tMB2 - MBase2 ) == MB2Dim + 1 );
   }
 #endif

 }  // end( CheckMB2 )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline void BMinQuad::CheckGS( void )
{
 #if( CHECK_DS & 1 )
  cIndex_Set tB2 = Base2;
  cIndex_Set tMB2 = MBase2;
  for( Index i = 0 ; i < SpaceDim ; i++ )
   if( GS[ i ] & kIsIn )         // an existing variable
    if( GS[ i ] & kInB2 )        // in Base2[]
     assert( i == *(tB2++) );
    else
     assert( i == *(tMB2++) );  // in MBase2[]
   else                         // a non-existing variable
    assert( ( *tB2 != i ) && ( *tMB2 != i ) );
 #endif

 }  // end( CheckGS )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

#define CMP( x1 , x2 ) ABS( x1 ) > eD * std::max( ABS( x2 ) , HpNum( 1 ) )

inline void BMinQuad::CheckRA( void )
{
 #if( CHECK_DS & 2 )
  for( Index i = 0 ; i < NxtBIdx ; i++ )
   if( IsThere( i ) ) {
    cHpNum xy = RealAlfa[ i ];
    cHpNum xx = xy - Alfa[ i ] - GiTLB( i , bounds , Base2 , B2Dim );

    if( CMP( xx , xy ) )
     #if( LOG_BMQ )
      *BMQLog << endl << "Test failed: RealAlfa[" << i << "] = " << xx;
     #else
      assert( false );
     #endif
    }
 #endif

 }  // end( CheckRA )

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

inline void BMinQuad::CheckDS( void )
{
 CheckB2();

 CheckMB2();

 CheckGS();

 CheckRA();
 }

/*--------------------------------------------------------------------------*/
/*---------------------- End File BMinQuad.C -------------------------------*/
/*--------------------------------------------------------------------------*/
