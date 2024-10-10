/*--------------------------------------------------------------------------*/
/*--------------------------- File ColorTV.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Definition and implementation of the class ColorTV, which is a stepsizes 
 * rule (SR) to be used within the subgradient method (SM) [see SubGrad.h].
 * The class conforms to the interface defined by the class Stepsize [see
 * Stepsize.h]. For more details we refer to:
 *
 * F. Barahona and R. Anbil <em>The Volume Algorithm: Producing Primal
 * Solutions with a Subgradient Method</em> Math. Prog. 87(3), 385--399,
 * 2000
 *
 * \author Antonio Frangioni \n
 *         Department of Informatics  \n
 *         University of Pisa \n
 *
 * \author Enrico Gorgone \n
 *         Department of Informatics \n
 *         Free University of Brussels \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*-------------------- DEFINITIONS & IMPLEMENTATION ------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ColorTV
 #define __ColorTV  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Stepsize.h"

/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS  ----------------------------------*/
/*--------------------------------------------------------------------------*/

#define LOG_STP 0

/* If LOG_STP > 0, the ColorTV class produces a log of its activities on the
   ostream object and at the "level of verbosity" set with the method
   SetSTPLog() [see below]. */

#if( LOG_STP )
 #define STPLOG( l , x ) if( STPLLvl > l ) *STPLog << x
 #define STPLOG2( l , c , x ) if( ( STPLLvl > l ) && c ) *STPLog << x
#else
 #define STPLOG( l , x )
 #define STPLOG2( l , c , x )
#endif

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** Definition of the class ColorTV. This class implements the target stepsize
    used in the original Volume algorithm. The method is based on classifying
    the iterations based on the obtained improvement of the function Fi(). */

class ColorTV : public Stepsize
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
    The ColorTV defines one main public types: */

 enum condition { green ,  ///< good step
		  yellow , ///< not bad, not good step
		  red      ///< bad step
                  };

/** @} ---------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** Constructor of the class. Since the constructor of ColorTV is executed
   after the one of Stepsize, the following parameters specific for ColorTV
   have to be found in the stream <em>after</em> those of the base class
   [see the comments to the constructor of Stepsize]:

     -# HpNum BetaZero      [1] initial value of beta

     -# Index greentestinvl [1] how many consecutive green iterations
                                are allowed before increasing beta

     -# Index yellowtestinvl [400] how many consecutive yellow iterations
                             are allowed before decreasing beta

     -# Index redtestinvl    [10] how many consecutive red iterations are
                             allowed before decreasing beta */

 ColorTV( SubGrad *slvr , std::istream *iStrm = 0 )
  : Stepsize( slvr , iStrm )
 {
  DfltdSfInpt( iStrm , BetaZero , HpNum( 1 ) );  // initial value of beta 
  DfltdSfInpt( iStrm , greentestinvl , Index( 1 ) );
  DfltdSfInpt( iStrm , yellowtestinvl , Index( 400 ) );
  DfltdSfInpt( iStrm , redtestinvl , Index( 10 ) );
  }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

 inline void SetSTPLog( ostream *outs = 0 , const char lvl = 0 );

/*--------------------------------------------------------------------------*/

 inline void Format( void ) ;

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR STEPSIZE COMPUTATION --------------------*/
/*--------------------------------------------------------------------------*/

 inline void NewStep( void );

/*--------------------------------------------------------------------------*/

 bool NeedsdkM1Gk( void ) { return( true ); }

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 inline bool UpdateTargetLevel( void );

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 inline condition coloring( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
/*--------------------------------------------------------------------------*/

 Index greentestinvl;   ///< maximum number of consecutive green iteration
 Index yellowtestinvl;  ///< maximum number of consecutive yellow iteration
 Index redtestinvl;     ///< maximum number of consecutive red iteration

 Index lastgreeniter;   ///< the last green iteration
 Index lastyellowiter;  ///< the last yellow iteration
 Index lastrediter;     ///< the last red iteration

 HpNum BetaZero;
 HpNum LwrBnd;     ///< lower bound
 HpNum lastvalue;  ///< the previous function value

 HpNum FiLmbd;     ///< current FiLambda
 Index NrIter;

/*--------------------------------------------------------------------------*/

 };  // end( class ColorTV )

/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

inline void ColorTV::Format( void )
{
 Stepsize::Format();

 lastgreeniter = lastyellowiter = lastrediter = 0;
 LwrBnd = - Inf< HpNum >();           // lower bound, target level are unknown
 lastvalue = FiLmbd = Inf< HpNum >(); // \f$ \bar{f}_{i-1} \f$ and \f$ f_i \f$
                                      // are unknown
 Beta = BetaZero;
 NrIter = 0;
 }

/*--------------------------------------------------------------------------*/

inline void ColorTV::NewStep( void )
{
 // get the function value - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 FiLmbd = ReadFkVal();

 if( FiLmbd == Inf< HpNum >() )
  throw( NDOException( "ColorTV::GetStepsize: this should not happen" ) );

 // eventually, improve the target level  - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 UpdateTargetLevel();

 // which is it the color?  - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 coloring();
 NrIter++;

 // print the level   - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 STPLOG( 1 , std::endl << "           " << " level  = " << -FiLev
	               << " ~  beta = " << Beta << std::endl << "           "
	 );

 }  // end( ColorTV::NewStep )

/*--------------------------------------------------------------------------*/

inline void ColorTV::SetSTPLog( ostream *outs , const char lvl )
{
 Stepsize::SetSTPLog( outs , lvl );

 #if( LOG_STP )
  if( STPLLvl > 1 ) {
   *STPLog << std::endl << "ColorTV: ~ BetaZero = " << BetaZero
	   << " ~ Color (g , y , r) = (" << greentestinvl << " , "
	   << yellowtestinvl << " , " << redtestinvl << ")" << std::endl;
   }
 #endif
 }

/*--------------------------------------------------------------------------*/

inline bool ColorTV::UpdateTargetLevel( void )
{
 bool TLchgd = false;

 // the target level is updated after a change of the lower bound  - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 cHpNum LwrBnd_ = GetOracle()->GetLowerBound();
 if( LwrBnd_ > LwrBnd )
  LwrBnd = LwrBnd_;

 // some exception - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( LwrBnd == -Inf< HpNum >() )
  throw( NDOException(
		    "ColorTV::UpdateTargetLevel: no lower bound is given" ) );

 // target level updating  - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( FiLev < -10 ) || ( FiLev > 10 ) ) {  // nicely away from 0
  if( FiLmbd <= FiLev + ABS( FiLev ) * 0.05 ) {
   if( ( FiLmbd < 10 ) && ( FiLmbd > - 9.5 ) )
    FiLev= - 10;
   else {
    FiLev -= 0.025 * ABS( FiLev  );
    FiLev = std::min( FiLev , FiLmbd - 0.05 * ABS( FiLmbd ) );
    }
   TLchgd = true;
   }
  }
 else                                       // near 0, must be careful
  if( FiLmbd < 10 ) {
   FiLev = ( FiLmbd < - 9.5 ) ? - 10 :
            ( FiLmbd > 0 ? FiLmbd * 0.95 : FiLmbd * 1.05 );
   TLchgd = true;
   }

 // failure of the target level updating - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( FiLev < LwrBnd ) {  // change the level
  FiLev = LwrBnd;
  TLchgd = true;
  }

 return( TLchgd );

 }  // end( ColorTV::UpdateTargetLevel )

/*--------------------------------------------------------------------------*/

inline ColorTV::condition ColorTV::coloring( void )
{
 Index ConColor;  // number of consecutive iteration with the same color

 // green, yellow, red - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 condition lastswing;
 if( ( -GetdGk() <= -1e-6 ) && ( lastvalue - FiLmbd >=
	       1e-6 * std::max( ABS( Solver->ReadBestFiVal() ) , 1.0 ) ) )  {
  lastswing = green;
  lastgreeniter = NrIter;
  }
 else
  if( ( -GetdGk() > -1e-6 ) && ( lastvalue > FiLmbd ) ) {
   lastswing = yellow;
   lastyellowiter = NrIter;
   }
  else {
   lastswing = red;
   lastrediter = NrIter;
   }

 // change beta, if needed - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 STPLOG( 1 , std::endl << "           " );

 switch( lastswing ) {
  case( green ):  // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ConColor = NrIter - std::max( lastyellowiter , lastrediter );
   STPLOG( 1 , " # green =  " << ConColor << " ~ ScPr = "
	       << GetdGk() <<  " ~ DeltaFi = " << lastvalue - FiLmbd );
   if( ( ConColor >= greentestinvl ) && ( Beta < 2 ) ) {
    lastgreeniter = lastyellowiter = lastrediter = NrIter;
    Beta *=  2.0;
    }
   break;

  case( yellow ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ConColor = NrIter - std::max( lastgreeniter , lastrediter );
   STPLOG( 2 , " # yellow =  " << ConColor << " ~ ScPr = "
	       << GetdGk() <<  " ~ DeltaFi = " << lastvalue - FiLmbd );
   if( ConColor >= yellowtestinvl ) {
    lastgreeniter = lastyellowiter = lastrediter = NrIter;
    Beta *= 1.1;
    }
   break;

  case( red ):  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   ConColor = NrIter - std::max( lastgreeniter , lastyellowiter );
   STPLOG( 1 , " # red =  " << ConColor << " ~ ScPr = " << GetdGk() );
   if( ( ConColor >= redtestinvl ) && ( Beta > HpNum( 5e-4 ) ) ) {
    lastgreeniter = lastyellowiter = lastrediter = NrIter;
    Beta *= 0.67;
    }
   break;

  }  // end( switch( lastswing ) )- - - - - - - - - - - - - - - - - - - - - -

 // update the last value - - - - - - - - - - - - - - - - - - - - - - - - - -

 lastvalue = Solver->ReadFiVal();
 return( lastswing );

 }  // end( ColorTV::coloring )

/*--------------------------------------------------------------------------*/

};  // end( namespace NDO_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ColorTV.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File ColorTV.h -----------------------------*/
/*--------------------------------------------------------------------------*/
