/*--------------------------------------------------------------------------*/
/*---------------------------- File MPTester.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition of the MPTester class, which implements the MPSolver interface.
 * The class MPTester is a special MPSolver that compares and logs the runs
 * of two MPSolver: Master and Slave. MPTester redirects every call to its
 * methods to the respective methods of Master and Slave, assuming that
 * Master is a "correct" MPSolver and Slave is an "uncertain" MPSolver. The
 * MPTester behavior is equivalent to Master behavior, but checks are done
 * to signal if Slave's behaviour deviates from Master's. 
 *
 * \author Luigi Poderico \n
 *         http:\\poderico.supereva.it \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MPTester
 #define __MPTester  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MPSolver.h"
#include <algorithm>

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

 using namespace std;

 using namespace OPTtypes_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CLASSES ------------------------------------*/
/*--------------------------------------------------------------------------*/
/** The MPTester class implements the MPSolver interface. An object of class
 * MPTester compares and logs the runs of two MPSolver: Master and Slave.
 * MPTester redirects every call to itt methods to the respective methods of
 * Master and Slave, assuming that Master is a "correct" MPSolver and Slave is
 * an "uncertain" MPSolver. The MPTester behavior is equivalent to Master
 * behavior, but checks are done to signal if Slave's behaviour deviates from
 * Naster's. */

class MPTester : public MPSolver
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/

 MPTester( MPSolver * master , MPSolver * slave )
  : fMster( master ) , fSlave( slave )
 {
  if( ! fMster )
   throw( std::invalid_argument( "null master MPSolver" ) );
  if( ! fSlave )
   throw( std::invalid_argument( "null slave MPSolver" ) );
  if( ! MPLog )
   MPLog = &std::cout;
  }

/*--------------------------------------------------------------------------*/

 MPSolver * get_master( void ) { return( fMster ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 MPSolver * get_slave( void ) { return( fSlave ); }

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

 virtual ~MPTester() { delete fMster; delete fSlave; };

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 void SetDim( cIndex MxBSz = 0 , FiOracle * Oracle = nullptr ,
	      const bool UsAvSt = false ) override final
 {
  fMster->SetDim( MxBSz , Oracle , UsAvSt );
  fSlave->SetDim( MxBSz , Oracle , UsAvSt );
  if( MxBSz )
   CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void Sett( cHpNum tt = 1 ) override final
 {
  fMster->Sett( tt );
  fSlave->Sett( tt );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void SetPar( const int wp , cHpNum value ) override final
 {
  fMster->SetPar( wp , value );
  fSlave->SetPar( wp , value );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void SetThreads( int nthreads ) override final
 {
  fMster->SetThreads( nthreads );
  fSlave->SetThreads( nthreads );
  }

/*--------------------------------------------------------------------------*/

 void SetLowerBound( cHpNum LwBnd = - Inf< HpNum >() ,
		     cIndex wFi = Inf< Index >() ) override final
 {
  fMster->SetLowerBound( LwBnd , wFi );
  fSlave->SetLowerBound( LwBnd , wFi );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void CheckIdentical( const bool Chk = true ) override final
 {
  fMster->CheckIdentical( Chk );
  fSlave->CheckIdentical( Chk );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void SetMPLog( std::ostream * outs = NULL ,
		const char lvl = 0 ) override final
 {
  fMster->SetMPLog( outs , lvl );
  fSlave->SetMPLog( outs , lvl );
  MPSolver::SetMPLog( outs , lvl );
  if( ! MPLog )
   MPLog = &std::cout;
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void SetMPTime( const bool TimeIt = true ) override final
 {
  fMster->SetMPTime( TimeIt );
  fSlave->SetMPTime( TimeIt );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

 MPStatus SolveMP( void ) override final
 {
  CheckProblem();

  MPStatus master = fMster->SolveMP();
  MPStatus slave = fSlave->SolveMP();

  if( master != slave )
   abort();

  return( master );
  }

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

 HpNum ReadFiBLambda( cIndex wFi = Inf< Index >() ) override final
 {
  // Check primal solution
  // Readd(true);

  HpNum master = fMster->ReadFiBLambda( wFi );
  HpNum slave = fSlave->ReadFiBLambda( wFi );

  if( wFi == Inf< Index >() )
   CompareDouble( master , slave , "ReadFiBLambda" );
  else
   if( wFi > fMster->GetNrFi() )
    CompareDouble( master , slave , "ReadFiBLambda( all - 0 )" );
   else {
    std::string mess( "ReadFiBLambda( " + std::to_string( wFi ) + " )" );
    CompareDouble( master , slave , mess.data() );
    }

  return( master );
  }

/*--------------------------------------------------------------------------*/

 HpNum ReadDt( cHpNum tt = 1 ) override final
 {
  HpNum master = fMster->ReadDt( tt );
  HpNum slave = fSlave->ReadDt( tt );

  CompareDouble( master , slave , "ReadDt" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 HpNum ReadSigma( cIndex wFi = Inf< Index >() ) override final
 {
  HpNum master = fMster->ReadSigma( wFi );
  HpNum slave = fSlave->ReadSigma( wFi );

  if( wFi > fMster->GetNrFi() )
   CompareDouble( master , slave , "ReadSigma" );
  else {
   std::string mess( "ReadSigma( " + std::to_string( wFi ) + " )" );
   CompareDouble( master , slave , mess.data() );
   }

  return( master );
  }

/*--------------------------------------------------------------------------*/

 HpNum ReadDStart( cHpNum tt = 1 ) override final
 {
  HpNum master = fMster->ReadDStart( tt );
  HpNum slave = fSlave->ReadDStart( tt );

  CompareDouble( master , slave , "ReadDStart" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 cLMRow Readd( bool Fulld = false ) override final
 {
  cLMRow master = fMster->Readd( Fulld );
  cLMRow slave = fSlave->Readd( Fulld );

  CompareInt( fMster->GetCrrSGLen() , fSlave->GetCrrSGLen() ,
		    "GetCrrSGLen" );
  CompareVector( master , slave , fMster->GetCrrSGLen() , "d" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 void ReadZ( LMRow tz , cIndex_Set &I , Index &D ,
	     cIndex wFi = Inf< Index >() ) override final
 {
  static std::vector< LMNum > slave;
  slave.resize( fSlave->GetCrrSGLen() );
  cIndex_Set slaveI;
  Index slaveD;

  fMster->ReadZ( tz , I , D , wFi );
  fSlave->ReadZ( slave.data() , slaveI , slaveD , wFi );

  CompareInt( fMster->GetCrrSGLen() , fSlave->GetCrrSGLen() ,
		    "GetCrrSGLen" );
  CompareInt( D , slaveD , "ReadZ-D" );
  CompareVector( tz , slave.data() , D , "z" );
  }

/*--------------------------------------------------------------------------*/

 cHpRow ReadMult( cIndex_Set &I , Index &D , cIndex wFi = Inf< Index >() ,
		  const bool IncldCnst = true ) override final
 {
  cHpRow master = fMster->ReadMult( I , D , wFi , IncldCnst );
		
  Index slaveD;
  cIndex_Set slaveI;
  cHpRow slave = fSlave->ReadMult( slaveI , slaveD , wFi , IncldCnst );

  CompareInt( fMster->GetMaxBSize() , fSlave->GetMaxBSize() ,
		    "GetMaxBSize" );

  if( fMster->GetMaxBSize()==0 )
   return( master );
		
  // check for mult vectors
  static std::vector< HpNum > denseD;
  static std::vector< HpNum > denseSlaveD;
  denseD.resize( fMster->GetMaxBSize() );
  denseSlaveD.resize( fSlave->GetMaxBSize() );
  
  size_t i;
  for( i = 0 ; i < fMster->GetMaxBSize() ; ++i )
   denseD[ i ] = denseSlaveD[ i ] = 0;

  for( i = 0 ; i < D ; ++i )
   denseD[ I[ i ] ] = master[ i ];

  for( i = 0 ; i < slaveD ; ++i )
   denseSlaveD[ slaveI[ i ] ] = slave[ i ];

  if( wFi > fMster->GetNrFi() )
  CompareVector( denseD.data() , denseSlaveD.data() ,
		       fMster->GetMaxBSize() , "Mult" );
  else {
   std::string mess( "Mult( " + std::to_string( wFi ) + " )" );
   CompareVector( denseD.data() , denseSlaveD.data() ,
			fMster->GetMaxBSize() , mess.data() );
   }

  return( master );
  }

/*--------------------------------------------------------------------------*/

 HpNum ReadLBMult( cIndex wFi = Inf< Index >() ) override final
 {
  HpNum master = fMster->ReadLBMult( wFi );
  HpNum slave = fSlave->ReadLBMult( wFi );

  if( wFi > fMster->GetNrFi() )
   CompareDouble( master , slave , "LBMult" );
  else {
   std::string mess( "LBMult( " + std::to_string( wFi ) + " )" );
   CompareDouble( master , slave , mess.data() );
   }

  return( master );
  }

/*--------------------------------------------------------------------------*/

 HpNum ReadGid( cIndex Nm = Inf< Index >() ) override final
 {
  HpNum master = fMster->ReadGid( Nm );
  HpNum slave = fSlave->ReadGid( Nm );

  CompareDouble( master , slave , "Gid" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 void MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau ) override final
 {
  fMster->MakeLambda1( Lmbd , Lmbd1 , Tau );

  static std::vector< HpNum > sLmbd1;
  sLmbd1.resize( fSlave->GetCrrSGLen() );
  CompareInt( fMster->GetCrrSGLen() , fSlave->GetCrrSGLen() ,
		    "GetCrrSGLen" );
  fSlave->MakeLambda1( Lmbd , sLmbd1.data() , Tau );

  CompareVector( Lmbd1 , sLmbd1.data() , fMster->GetCrrSGLen() ,
		       "Lambda1" );
  }

/*--------------------------------------------------------------------------*/

 void SensitAnals( HpNum &lp , HpNum &cp ) override final
 {
  fMster->SensitAnals( lp , cp );
		
  HpNum slp, scp;
  fSlave->SensitAnals( slp , scp );

  CompareDouble( lp , slp , "SensitAnals-lp" );
  CompareDouble( cp , scp , "SensitAnals-cp" );
  }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/

 Index GetMaxBSize( void )
 {
  Index master = fMster->GetMaxBSize();
  Index slave = fSlave->GetMaxBSize();

  CompareInt( master , slave , "MaxBSize" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 Index GetMaxSGLen( void )
 {
  Index master = fMster->GetMaxSGLen();
  Index slave = fSlave->GetMaxSGLen();

  CompareInt( master , slave , "MaxSGLen" );
  return( master );
  }

 /*--------------------------------------------------------------------------*/

 Index GetCrrSGLen( void )
 {
  Index master = fMster->GetCrrSGLen();
  Index slave = fSlave->GetCrrSGLen();

  CompareInt( master , slave , "CrrSGLen" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 Index GetNrFi( void )
 {
  Index master = fMster->GetNrFi();
  Index slave = fSlave->GetNrFi();

  CompareInt( master , slave , "NrFi" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 Index BSize( cIndex wFi = Inf< Index >() ) override final
 {
  Index master = fMster->BSize( wFi );
  Index slave = fSlave->BSize( wFi );

  CompareInt( master , slave , "BSize" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 Index BCSize( cIndex wFi = Inf< Index >() ) override final
 {
  Index master = fMster->BCSize( wFi );
  Index slave = fSlave->BCSize( wFi );

  CompareInt( master , slave , "BCsize" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 Index MaxName( cIndex wFi = Inf< Index >() ) override final
 {
  Index master = fMster->MaxName( wFi );
  Index slave = fSlave->MaxName( wFi );

  CompareInt( master , slave , "MaxName" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 Index WComponent( cIndex i ) override final
 {
  Index master = fMster->WComponent( i );
  Index slave = fSlave->WComponent( i );

  CompareInt( master , slave , "WComponent" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 bool IsSubG( cIndex i ) override final
 {
  bool master = fMster->IsSubG( i );
  bool slave = fSlave->IsSubG( i );

  CompareBool( master , slave , "IsSubG" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 Index NumNNVars( void ) override final
 {
  Index master = fMster->NumNNVars();
  Index slave = fSlave->NumNNVars();

  CompareInt( master , slave , "NumNNVars" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 Index NumBxdVars( void ) override final
 {
  Index master = fMster->NumBxdVars();
  Index slave = fSlave->NumBxdVars();

  CompareInt( master , slave , "NumBxdVars" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 bool IsNN( cIndex i ) override final
 {
  bool master = fMster->IsNN(i);
  bool slave = fSlave->IsNN(i);

  CompareBool( master , slave , "IsNN" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 cHpRow ReadLinErr( void ) override final
 {
  cHpRow master = fMster->ReadLinErr();
  cHpRow slave = fSlave->ReadLinErr();

  bool eq = true;
  cIndex mBS = fMster->GetMaxBSize();
  CompareInt( mBS , fSlave->GetMaxBSize() , "GetMaxBSize" );
  for( int i = 0 ; i < mBS ; ++i )
   if( fMster->WComponent( i ) < Inf< Index >() )
    if( ! CompareDouble( master[ i ] , slave[ i ] ) ) {
     eq = false;

     *MPLog << std::endl << "LinErr[ " << i << " ] = "
	    << master[ i ] << ", " << slave[ i ];
     }
     
  return( master );
  }

/*--------------------------------------------------------------------------*/

 HpNum ReadLowerBound( cIndex wFi = Inf< Index >() ) override final
 {
  HpNum master = fMster->ReadLowerBound( wFi );
  HpNum slave = fSlave->ReadLowerBound( wFi );

  if( wFi > fMster->GetNrFi() )
   CompareDouble( master , slave , "GlobalLowerBound" );
  else {
   std::string mess( "LowerBound( " + std::to_string( wFi ) + " )" );
   CompareDouble( master , slave , mess.data() );
   }

  return( master );
  }

/*--------------------------------------------------------------------------*/

 HpNum EpsilonD( void ) override final
 {
  HpNum master = fMster->EpsilonD();
  HpNum slave = fSlave->EpsilonD();

  CompareDouble( master , slave , "EpsilonD" );
  return( master );
  }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/

 SgRow GetItem( cIndex wFi = Inf< Index >() ) override final
 {
  fMsterItem = fMster->GetItem( wFi );
  fSlaveItem = fSlave->GetItem( wFi );

  CompareInt( fMster->GetCrrSGLen() , fSlave->GetCrrSGLen() ,
		    "CrrSGLen" );
  return( fMsterItem );
  }

/*--------------------------------------------------------------------------*/

 void SetItemBse( cIndex_Set SGBse = 0 , cIndex SGBDm = 0 ) override final
 {
  memcpy( fSlaveItem , fMsterItem , sizeof( SgNum ) * SGBDm );
		
  fMster->SetItemBse( SGBse , SGBDm );
  fSlave->SetItemBse( SGBse , SGBDm );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 Index CheckSubG( cHpNum DFi , cHpNum Tau , HpNum &Ai , HpNum &ScPri )
  override final
 {
  HpNum slaveAi = Ai;
  HpNum slaveScPri = ScPri;
		
  Index master = fMster->CheckSubG( DFi , Tau , Ai , ScPri );
  Index slave = fSlave->CheckSubG( DFi , Tau , slaveAi , slaveScPri );

  CompareInt( master , slave , "CheckSubG ret" );
  CompareDouble( Ai , slaveAi , "CheckSubG Ai" );
  CompareDouble( ScPri , slaveScPri , "CheckSubG ScPri" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 Index CheckCnst( HpNum &Ai , HpNum &ScPri , cHpRow CrrPnt ) override final
 {
  HpNum slaveAi = Ai;
  HpNum slaveScPri = ScPri;
		
  Index master = fMster->CheckCnst( Ai , ScPri , CrrPnt );
  Index slave = fSlave->CheckCnst( slaveAi , slaveScPri , CrrPnt );

  CompareInt( master , slave , "CheckCnst ret" );
  CompareDouble( Ai , slaveAi , "CheckCns Ai" );
  CompareDouble( ScPri , slaveScPri , "CheckCnst ScPri" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 bool ChangesMPSol( void ) override final
 {
  bool master = fMster->ChangesMPSol();
  bool slave = fSlave->ChangesMPSol();

  CompareBool( master , slave , "ChangesMPSol" );
  return( master );
  }

/*--------------------------------------------------------------------------*/

 void SetItem( cIndex Nm = Inf< Index >() ) override final
 {
  fMster->SetItem( Nm );
  fSlave->SetItem( Nm );
  }

/*--------------------------------------------------------------------------*/

 void SubstItem( cIndex Nm ) override final
 {
  fMster->SubstItem( Nm );
  fSlave->SubstItem( Nm );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void RmvItem( cIndex i ) override final
 {
  fMster->RmvItem( i );
  fSlave->RmvItem( i );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void RmvItems( void ) override final
 {
  fMster->RmvItems();
  fSlave->RmvItems();
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void SetActvSt( cIndex_Set AVrs = 0 , cIndex AVDm = 0 ) override final
 {
  fMster->SetActvSt( AVrs , AVDm );
  fSlave->SetActvSt( AVrs , AVDm );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void AddActvSt( cIndex_Set Addd , cIndex AdDm , cIndex_Set AVrs )
  override final
 {
  fMster->AddActvSt( Addd , AdDm , AVrs );
  fSlave->AddActvSt( Addd , AdDm , AVrs );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void RmvActvSt( cIndex_Set Rmvd , cIndex RmDm , cIndex_Set AVrs )
  override final
 {
  fMster->RmvActvSt( Rmvd , RmDm , AVrs );
  fSlave->RmvActvSt( Rmvd , RmDm , AVrs );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void AddVars( cIndex NNwVrs ) override final
 {
  fMster->AddVars( NNwVrs );
  fSlave->AddVars( NNwVrs );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void RmvVars( cIndex_Set whch = 0 , Index hwmny = 0 ) override final
 {
  fMster->RmvVars( whch , hwmny );
  fSlave->RmvVars( whch , hwmny );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void ChgAlfa( cHpRow DeltaAlfa ) override final
 {
  fMster->ChgAlfa( DeltaAlfa );
  fSlave->ChgAlfa( DeltaAlfa );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void ChgAlfa( cHpRow NewAlfa , cIndex wFi ) override final
 {
  fMster->ChgAlfa( NewAlfa , wFi );
  fSlave->ChgAlfa( NewAlfa , wFi );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void ChgAlfa( cIndex i , cHpNum Ai ) override final
 {
  fMster->ChgAlfa( i , Ai );
  fSlave->ChgAlfa( i , Ai );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void ChangeCurrPoint( cLMRow DLambda , cHpRow DFi ) override final
 {
  fMster->ChangeCurrPoint( DLambda , DFi );
  fSlave->ChangeCurrPoint( DLambda , DFi );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void ChangeCurrPoint( cHpNum Tau , cHpRow DFi ) override final
 {
  fMster->ChangeCurrPoint( Tau , DFi );
  fSlave->ChangeCurrPoint( Tau , DFi );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/

 void ChgSubG( cIndex strt , Index stp , cIndex wFi ) override final
 {
  fMster->ChgSubG( strt , stp , wFi );
  fSlave->ChgSubG( strt , stp , wFi );
  CheckProblem();
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 bool CompareDouble( double master , double slave )
 {
  if( master == Inf< HpNum >() )
   return( slave == Inf< HpNum >() );

  if( master == - Inf< HpNum >() )
   return( slave == - Inf< HpNum >() );

  if( ( slave == Inf< HpNum >() ) || ( slave == - Inf< HpNum >() ) )
   return( false );

  const double kMinPrecision = 1e-4;

  return( std::abs( master - slave ) <= kMinPrecision *
	  std::max( std::abs( master ) , std::max( double( 1 ) ,
						   std::abs( slave ) ) ) );
  }

/*--------------------------------------------------------------------------*/

 bool CompareDouble( double master , double slave ,
		     const char * const caller )
 {
  if( CompareDouble( master , slave ) )
   return( true );

  if( caller )
   *MPLog << std::endl << caller << ": " << master << ", " << slave;
  return( false );
  }

/*--------------------------------------------------------------------------*/

 bool CompareInt( unsigned int master , unsigned int slave ,
		  const char * const caller = 0 )
 {
  if( master == slave )
   return( true );

  if( caller )
   *MPLog << std::endl << caller << ": " << master << ", " << slave;
  return( false );
  }

/*--------------------------------------------------------------------------*/

 bool CompareBool( bool master , bool slave , const char * const caller )
 {
  if( master == slave )
   return( true );

  if( caller )
   *MPLog << std::endl << caller << ": " << master << ", " << slave;
  return( false );
  }

/*--------------------------------------------------------------------------*/

 bool CompareVector( const double * master , const double * slave , int len ,
		     const char * const caller )
 {
  bool eq = true;

  for( int i = 0 ; i < len ; ++i )
   if( ! CompareDouble( master[ i ] , slave[ i ] ) ) {
    eq = false;

    *MPLog << std::endl << *caller << "[ " << i << " ] = "
	   << master[ i ] << ", " << slave[ i ];
    }

  return( eq );
  }

/*--------------------------------------------------------------------------*/

 void CheckProblem( void )
 {
  // Check some simple conditions
  assert( NumBxdVars() >= NumNNVars() );
  Index i = GetCrrSGLen();
  while( i-- )
   IsNN( i );

  ReadLinErr();
  }

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 MPSolver *fMster;
 MPSolver *fSlave;

 SgRow fMsterItem;
 SgRow fSlaveItem;

/*--------------------------------------------------------------------------*/

 };  // end( class MPTester )

/*--------------------------------------------------------------------------*/

 }  //end( namespace( NDO_di_unipi_it ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* MPTester.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File MPTester.h -----------------------------*/
/*--------------------------------------------------------------------------*/
