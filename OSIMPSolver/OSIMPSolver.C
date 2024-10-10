/*--------------------------------------------------------------------------*/
/*--------------------------- File OSIMPSolver.C ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the OSIMPSolver class, which solves Master Problems for
 * Bundle algorithms using a generic OSISolverInterface object. This class
 * conforms to the interface defined by the class MPSolver [see MPSolver.h].
 *
 * \author Antonio Frangioni (original idea & implementation)\n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Enrico Gorgone (implementation)\n
 *         Dipartimento di Informatica \n
 *         Universita' della Calabria \n
 *
 * \author Andrea Nerli (implementation)\n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------------- MACROS ----------------------------------*/
/*--------------------------------------------------------------------------*/
/* Unfortunately, OsiSolverInterface objects do not do a bunch of things that
 * underlying solvers actually do, among which
 *
 * - dealing with quadratic problems (hence, using a quadratic stabilizing
 *   term);
 *
 * - stopping after a given time.
 *
 * This macro allows to specify that the OsiSolverInterface actually is of
 * some specific type, so that methods specific of the interface of the
 * solver can be used to do these. It's a bad hack (it often is when
 * dynamic_cast<>() and static_cast<>() are involved), but it's the least
 * worse we can do while allowing to actually use these features: it's
 * OsiSolverInterface's fault, not ours.
 *
 * In particular, since we know exactly the type of the underlying
 * OsiXXXSolverInterface when we cast, it would seem that we could
 * static_cast<>() instead of dynamic_cast<>() them; yet, this is not
 * possible since OsiSolverInterface is a virtual base class for
 * the interested OsiXXXSolverInterface (for which reason is unclear to me).
 *
 * Note that, at least, OSIMPSolver.h does not depend on the switch.
 *
 * Possible values are:
 *
 * 0  ==> any OsiSolverInterface, so none of these features can be used;
 *
 * 1  ==> OsiCpxSolverInterface;
 *
 * 2  ==> OsiGrbSolverInterface. */

#ifndef WHICH_OSI_MP
 #define WHICH_OSI_MP 2
#endif

/*--------------------------------------------------------------------------*/
/* If OSIMPSOLVERLOG > 0, the OSIMPSolver class produces a log of its
 * activities on the ostream object and at the "level of verbosity" set with
 * the method SetMPLog(). */

#define OSIMPSOLVERLOG 0

#if( OSIMPSOLVERLOG )
 #define MSG( l , m ) if( MPLLvl > l ) (*MPLog) << m << flush
#else
 #define MSG( l , m )
#endif

/*--------------------------------------------------------------------------*/
/* If CHECK_DS > 0, various data structures are checked for correctness
 * during the run of the algorithm, tremendously slowing down the algorithm
 * but allowing to debug the thing.
 * What data structures are checked is coded bit-wise in CHECK_DS:
 *
 * bit 0 (+ 1)  =>  the various dict_* are tested */

#define CHECK_DS 0

/*--------------------------------------------------------------------------*/
// true if i was in the last Master Problem
#define WasInMP( i ) ( ! ( wcomp[ i ] & LLBIndex ) )

// the component of i, assuming wcomp[ i ] < InINF
#define WComp( i ) ( wcomp[ i ] & ~IMask )

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"

#include "OSIMPSolver.h"

#include "OPTvect.h"

#include "NDOSlver.h"

#include <algorithm>
#include <numeric>

#if WHICH_OSI_MP == 1
 #include "OsiCpxSolverInterface.hpp"
 #include "ilcplex/cplex.h"
#elif WHICH_OSI_MP == 2
 #include "OsiGrbSolverInterface.hpp"
 #include "gurobi_c++.h"
#endif

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------------- AUXILIARY TYPES -----------------------------*/
/*--------------------------------------------------------------------------*/

class DerivedHandler : public CoinMessageHandler {
 public:
  virtual int print( void ) { return( 0 ); }
 };

/*--------------------------------------------------------------------------*/
/*------------------------------- CONSTANTS --------------------------------*/
/*--------------------------------------------------------------------------*/

static cIndex InINF = Inf< Index >();
static cLMNum LMINF = Inf< LMNum >();
static cHpNum HpINF = Inf< HpNum >();

// an Index with the most significant bit only == 1
static cIndex LBIndex = 1 << ( numeric_limits< Index >::digits - 1 );
// an Index with the second-most significant bit only == 1
static cIndex LLBIndex = LBIndex >> 1;
// masking the two most significant bits in an Index
static cIndex IMask = LBIndex + LLBIndex;

static double defOptEps = 1e-10;
static double defFsbEps = 1e-10;

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTRUCTOR -------------------------------*/
/*--------------------------------------------------------------------------*/

OSIMPSolver::OSIMPSolver( std::istream * iStrm ) : MPSolver()
{
 // reset the pointers to OsiSolver and FiOracle object  - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 osiSlvr = nullptr;
 FIO = nullptr;

 // set the default algorithm parameters - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 useactiveset = false;  // no active set
 first = true;          // no problem has been solved so far
 checkID = false;       // by default don't check for identical items

 stab = unset;	        // the stabilization term is unknown
 t = 1;			// the proximity parameter is set to a standard value
 MaxTime = 0;           // no time limit
 OptEps = defOptEps;    // default optimality parameter
 FsbEps = defFsbEps;    // default feasibility parameter
    
 NRCall = 0;            // number of calls to SolveMP is zero
    
 // the master is empty- - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MaxBSize = 0;		// max bundle is zero
 NSubG = NConst = 0;    // bundle is empty

 dict_item = nullptr;   // no items are in the master problem
 item_maxname = 0;
 dict_slack = nullptr;
 dict_stab = nullptr;

 wcomp = nullptr;
 weasy = nullptr;
 NrEasy = 1;

 NewItem = 0;	        // new item information is empty
 NewItemBse = 0;
 NewItemBDm = 0;

 NNVars = BxdVars = 0;
 Upper = Lower = nullptr;
 LwrBnds = nullptr;

 tempHP = nullptr;       // temporary buffers of HpNum
 tempI = nullptr;        // and Index
 tempHP_size = 0;        // initially they are empty
 tempI_size = 0;

 Aset = nullptr;         // no active sets
 Asetdim = 0;

 algorithm = kAuto;       // automatic solution by default

 // no solutions is given  - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 comp_row = nullptr;
 comp_col = nullptr;

 RhoCol = nullptr;
 RhoColBse = nullptr;
 RhoColBDm = RhoColSgPos = 0;

 GiPerd = nullptr;

 #if( PRESERVE_OSI_SOLS )
  csol = rsol = rcst = nullptr;
  csols = rsols = 0;
 #endif

 derhand = new DerivedHandler;

 }  // end( OSIMPSolver::OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetDim( cIndex MxBSz , FiOracle *Oracle ,
			  const bool UsAvSt )
{
 MSG( 0 , "OSIMPSolver::SetDim()\n" );

 if( ! MxBSz ) {  // deallocate all its memory and quietly wait
  cleanup();      // for new instructions- - - - - - - - - - - - - - - - - - -
  return;
  }

 if( Oracle ) {  // allocate the memory for solving the problem- - - - - - - -
  if( ! osiSlvr )
   throw( NDOException( "OSIMPSolver::SetDim: osiSlvr must be set" ) );

  // discards all the previous settings and deallocate all the memory- - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cleanup();

  // set the problem as a maximum problem  - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->setObjSense( -1 );

  // ask FiOracle the instance dimension - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  MaxBSize = MxBSz;
  FIO = Oracle;
  MaxSGLen = Oracle->GetMaxNumVar();
  CrrSGLen = Oracle->GetNumVar();
  NrFi = Oracle->GetNrFi();
  MaxNZ = Oracle->GetMaxNZ();
  useactiveset = UsAvSt;

  // allocate the memory - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  weasy = new bool[ NrFi + 1 ];
  weasy[ 0 ] = true;  // the 0-th is very easy
  NrEasy = 1;

  NewItem = new SgNum[ MaxNZ ];

  NSubG = new Index[ NrFi + 1 ];
  VectAssign( NSubG , Index( 0 ) , NrFi + 1 );
  NConst = new Index[ NrFi + 1 ];
  VectAssign( NConst , Index( 0 ) , NrFi + 1 );

  comp_row = new Index[ NrFi + 1 ];
  comp_col = new Index[ NrFi + 1 ];
  dict_item = new Index[ MaxBSize ];
  VectAssign( dict_item , InINF , MaxBSize );
  dict_slack = new Index[ MaxSGLen ];
  VectAssign( dict_slack , InINF , MaxSGLen );
  wcomp = new Index[ MaxBSize ];
  VectAssign( wcomp , InINF , MaxBSize );

  Upper = new LMNum[ MaxSGLen ];
  Lower = new LMNum[ MaxSGLen ];

  item_maxname = 0;

  // create all the rows for static constraints and (partial) set of rows
  // for the dynamic constraints. The first constraint is relative to rho.
  // Note that variable r has been deleted since it's redundant. So the first
  // constraint is 1 - rho >= 0. All in all in the absence of total lower
  // bound case this constraint should be rho = 1 (the default case).
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->addRow( 0 , 0 , 0 , 1 , 1 );

  // add the component-specific rows - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // while each "hard" component only produces a single convexity constraint,
  // for each easy component one may have to generate up to 2 * (BNR + BNC)
  // constraints. the issue is that a component of the general form
  //
  //    max   ( c[ i ] - Lambda * A[ i ] ) x[ i ]
  //          r[ i ] <= B[ i ] x[ i ] <= p[ i ]
  //          l[ i ] <= x[ i ] <= u[ i ]
  //
  // in fact has to be put in the master problem as
  //
  //    max   ( c[ i ] - Lambda * A[ i ] ) x[ i ]
  //          \rho r[ i ] <= B[ i ] x[ i ] <= \rho p[ i ]
  //          \rho l[ i ] <= x[ i ] <= \rho u[ i ]
  //
  // with \rho the variable corresponding to the global lower bound. this
  // means that:
  //
  // - a "true" ranged constraint r[ i ] <= B[ i ] x[ i ] <= p[ i ],
  //   that is originally a *single* row, has to become *two* rows;
  //   however, this is only true if r[ i ] < p[ i ] and they are *both 
  //   nonzero and finite*. in fact, if r[ i ] == p[ i ] then the
  //   constraint is
  //  
  //          \B[ i ] x[ i ] - \rho p[ i ] == 0
  //
  //   i.e., a single equality (ranged constraint with 0 range). similarly,
  //   if (say) r[ i ] == 0  then the constraint is
  //
  //          0 <= \B[ i ] x[ i ] - \rho p[ i ] <= 0
  //
  //   i.e., again a single equality (ranged constraint with 0 range), while
  //   if (say) r[ i ] == -INF  then the constraint is
  //
  //          (-INF <=) \B[ i ] x[ i ] - \rho p[ i ] <= 0
  //
  //   i.e., a single inequality
  //
  // - simple bounds l[ i ] <= x[ i ] <= u[ i ] become *two* rows: however,
  //   again this is only true l[ i ] and u[ i ] are *both nonzero and
  //   finite* (the case l[ i ] == u[ i ] is not considered). in fact, if
  //   (say) 0 == l[ i ] < u[ i ]  then the constraint is
  //
  //          0 <= x[ i ] <= \rho u[ i ]
  //
  //   i.e., a single new equality is added while keeping the lower bound.
  //   if (say) 0 == l[ i ] < u[ i ] == +INF then the constraint is
  //
  //          0 <= x[ i ] (<= INF)
  //
  //   and hence no new inequalities are added
  //
  // the organization of the rows corresponding to component i is the
  // following:
  //
  // - they start at position comp_row[ i - 1 ]
  //
  // - the first rows are those "replacing" the bound constraints, between 0
  //   and 2 * BNC; they are in the same order as the variables (columns), and
  //   if a variable has two rows they are consecutive (first the lb one, then
  //   the ub one)
  //
  // - the following rows are described by the (temporary) dictionary rowpos:
  //   rowpos[ i - 1 ] is a vector as long as there are rows in the original
  //   B[ i ] matrix, and rowpos[ i - 1 ][ j ] is the index of the (first, if
  //   the row has been duplicated) corresponding row in the master problem
  //
  // - if an original row j produces two ones they are consecutive (first the
  //   lhs one, then the the rhs one) in the master problem
  
  // the temporary dictionary for the position in the master problem of (the
  // first of the two copies of) each row is constructed now
  std::vector< std::vector< Index > > rowpos( NrFi );
  
  for( Index i = 0 ; i++ < NrFi ; ) {
   comp_row[ i - 1 ] = osiSlvr->getNumRows();
   if( auto BNC = FIO->GetBNC( i ) ) {  // an easy component - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    weasy[ i ] = true;
    ++NrEasy;

    auto lbd = new double[ BNC ];
    auto ubd = new double[ BNC ];
    double * lhs = nullptr;
    double * rhs = nullptr;

    auto BNR = FIO->GetBNR( i );  // number of rows of matrix B[ i ]
    if( BNR ) {
     rowpos[ i - 1 ].resize( BNR );
     lhs = new double[ BNR ];
     rhs = new double[ BNR ];
     }

    // only ask for lhs and rhs if there are rows at all
    FIO->GetBDesc( i , 0 , 0 , 0 , lhs , rhs , 0 , lbd , ubd );

    // first create the rows replacing the variable bounds- - - - - - - - - -
    // ... but only if the corresponding bound is finite and nonzero
    //
    // l[ i ] <= x[ i ] becomes  0 <= x[ i ] - \rho l[ i ] (< INF)
    // x[ i ] <= u[ i ] becomes  (-INF) <= x[ i ] - \rho u[ i ] <= 0

    for( Index j = 0 ; j < BNC ; ++j ) {
     if( lbd[ j ] && ( lbd[ j ] > -Inf< double >() ) )
      osiSlvr->addRow( 0 , 0 , 0 , 0 , osiSlvr->getInfinity() );

     if( ubd[ j ] && ( ubd[ j ] < Inf< double >() ) )
      osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() , 0 );
     }

    delete[] ubd;
    delete[] lbd;

    if( ! BNR )  // if there are no rows
     continue;   // all done

    // then create the rows representing the original ones- - - - - - - - - -
    // for each original row for which both bounds are finite, nonzero and
    // different two ones are created, otherwise only one is
    //
    // r[ i ] <= B[ i ] x[ i ] becomes 0 <= B[ i ] x[ i ] - \rho r[ i ] < INF
    // B[ i ] x[ i ] <= p[ i ] becomes -INF < B[ i ] x[ i ] - \rho p[ i ] <= 0

    for( Index j = 0 ; j < BNR ; ++j ) {
     if( ( lhs[ j ] == -Inf< double >() ) && ( rhs[ j ] == Inf< double >() ) )
      throw( std::invalid_argument( "-INF <= row <= +INF not allowed" ) );

     rowpos[ i - 1 ][ j ] = osiSlvr->getNumRows();  // record position

     if( lhs[ j ] && ( lhs[ j ] > -Inf< double >() ) &&
	 rhs[ j ] && ( rhs[ j ] < Inf< double >() ) &&
	 ( lhs[ j ] < rhs[ j ] ) ) {  // two rows are created
      osiSlvr->addRow( 0 , 0 , 0 , 0 , osiSlvr->getInfinity() );
      osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() , 0 );
      }
     else {                           // one row is created
      double lb = lhs[ j ] > -Inf< double >() ? 0 : -Inf< double >();
      double ub = rhs[ j ] < Inf< double >() ? 0 : Inf< double >();
      osiSlvr->addRow( 0 , 0 , 0 , lb , ub );
      }
     }  // end( for( j ) )

    delete[] rhs;
    delete[] lhs;
    }
   else {  // a "hard" component - - - - - - - - - - - - - - - - - - - - - - -
    weasy[ i ] = false;
    osiSlvr->addRow( 0 , 0 , 0 , 0 , 0 );  // create the simplex constraint
    }
   }  // end( for( each component ) )

  // mark the beginning of subgradient-related rows- - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  comp_row[ NrFi ] = osiSlvr->getNumRows();
  
  // generate the subgradient-related rows - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( useactiveset ) { // all the constraints are inactive so far- - - - - - -
   MSG( 0 , "The dynamic inactive constraints have been created starting "
	    "from the row " << comp_row[ NrFi ] << std::endl );

   for( Index j = CrrSGLen ; j > 0 ; j-- )
    osiSlvr->addRow( 0 , 0 , 0 , - osiSlvr->getInfinity() ,
		                   osiSlvr->getInfinity() );
   }
  else {     // all the constraints are active (now and for ever)- - - - - - -
   MSG( 0 , "Dynamic constraints have been created starting from the row "
	     << comp_row[ NrFi ] << std::endl );
   for( Index j = CrrSGLen ; j > 0 ; j-- )
    osiSlvr->addRow( 0 , 0 , 0 , 0 , 0 );
   }

  // allocate and initialize dictionary and values of the column relative to
  // the special variable \rho. if a global lower bound exists, the
  // coefficient of \rho in the objective function is *-* the bound;
  // otherwise the coefficient is 0 and \rho is to fixed/ to 1. the
  // coefficients in the column of \rho are the RHS/LHS of the constraints
  // defining the easy components (comprised the replacements of finite and
  // nonzero bounds), which are added next, and the gradient of 0-th
  // component, that is, the b vecto, which will be provided later
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  RhoCol = new double[ comp_row[ NrFi ] + MaxSGLen ];
  RhoColBse = new int[ comp_row[ NrFi ] + MaxSGLen ];

  RhoColBse[ 0 ] = 0;
  RhoCol[ 0 ] = 1;
  RhoColBDm = 1;

  // allocate the memory for the temporary vectors tempI and tempHp
  resizeI( osiSlvr->getNumRows() );
  resizeHP( osiSlvr->getNumRows() );

  // add the component-specific columns- - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // each component i has some "special" variables (columns) at positions
  // comp_col[ i ] <= j < comp_col[ i + 1 ]. for difficult components this is
  // just the \gamma^i variable associated with the individual lower bound (if
  // any), which also doubles as "artificial slack for the phase 0", while for
  // "easy" components they are all the variables x[ i ]
  //
  // HOWEVER, NOTE THAT "comp_col[ i ] <= j < comp_col[ i + 1 ]" IS ONLY VALID
  // HERE WHEN THE DATA IS BEING CONSTRUCTED, BUT NOT IN GENERAL; this is
  // because in ChgSubG() the order of variables may be changed (due of dire
  // limitations about how OsiSolverInterface allows you to change the matrix
  // of the problem). Hence, while it will remain true that the variables of
  // component i start at comp_col[ i ], comp_col[ i ] > comp_col[ i + 1 ]
  // may happen (and as a result, there is no easy way to know how many
  // columns are associated to each easy component)
  
  for( Index i = 0 ; ++i <= NrFi ; ) {
   // mark the start of the variables of the component, be it easy or not
   comp_col[ i ] = osiSlvr->getNumCols();

   if( auto BNC = FIO->GetBNC( i ) ) {  // number of columns of matrix B[ i ]
                                        // - - - - - - - - - - - - - - - - - -
    // ... if nonzero this is an easy component; note that FiOracle returns
    // all the information in a sparse format, so the columns of the matrices
    //  A and B are in sparse format

    auto BNR = FIO->GetBNR( i );  // number of rows of matrix B[ i ]
    auto BNZ = FIO->GetBNZ( i );  // number of non-zeroes of matrix B[ i ]
    auto ANZ = FIO->GetANZ( i );  // number of non-zeroes of matrix A[ i ]

    // allocate the memory for the costs - - - - - - - - - - - - - - - - - - -

    auto cst = new double[ BNC ];

    // allocate the memory for the matrix B[ i ], the lhs/rhs, the ub/lb - - -

    auto lbd = new double[ BNC ];
    auto ubd = new double[ BNC ];

    int * Bbeg = nullptr;
    int * Bind = nullptr;
    double * Bval = nullptr;
    double * lhs = nullptr;
    double * rhs = nullptr;
    if( BNR ) {
     Bbeg = new int[ BNC + 1 ];
     Bind = new int[ BNZ ];
     Bval = new double[ BNZ ];
     lhs = new double[ BNR ];
     rhs = new double[ BNR ];
     }

    // ask FiOracle for the matrix B[ i ] and all the other stuff- - - - - - -

    FIO->GetBDesc( i , Bbeg , Bind , Bval , lhs , rhs , cst , lbd , ubd );

    // allocate the memory to describe the matrix A[ i ] - - - - - - - - - - -

    auto Abeg = new int[ BNC + 1 ];
    auto Aind = new int[ ANZ ];
    auto Aval = new double[ ANZ ];

    // ask FiOracle for the matrix A[ i ]: note that the FiOracle may give
    // partial information of the rows of matrix A[ i ], and precisely the
    // first CrrSGLen rows only- - - - - - - - - - - - - - - - - - - - - - - -

    FIO->GetADesc( i , Abeg , Aind , Aval );

    // write the columns relative to the easy component- - - - - - - - - - - -
    // meanwhile, fill the corresponding nonzero entries of RhoCol for the
    // bounds replacements only

    Index boundrow = comp_row[ i - 1 ];   // first row of bounds replacement

    for( Index j = 0 ; j < BNC ; j++ ) {  // for each variable j
     Index count = 0;                     // number of nonzeros in the column

     // static part I: box constraints - - - - - - - - - - - - - - - - - - - -

     double lbj;
     if( lbd[ j ] && ( lbd[ j ] > -Inf< double >() ) ) {
      // l[ i ] <= x[ i ] becomes  0 <= x[ i ] - \rho l[ i ] (< INF)
      tempI[ count ] = boundrow;
      tempHP[ count++ ] = 1;
      RhoColBse[ RhoColBDm ] = boundrow++;
      RhoCol[ RhoColBDm++ ] = - lbd[ j ];
      lbj = - osiSlvr->getInfinity();
      }
     else  // the bound remains a bound
      lbj = lbd[ j ] == 0 ? 0 : - osiSlvr->getInfinity();

     double ubj;
     if( ubd[ j ] && ( ubd[ j ] < Inf< double >() ) ) {
      // x[ i ] <= u[ i ] becomes  (-INF) <= x[ i ] - \rho u[ i ] <= 0
      tempI[ count ] = boundrow;
      tempHP[ count++ ] = 1;
      RhoColBse[ RhoColBDm ] = boundrow++;
      RhoCol[ RhoColBDm++ ] = - ubd[ j ];
      ubj = osiSlvr->getInfinity();
      }
     else  // the bound remains a bound
      ubj = ubd[ j ] == 0 ? 0 : osiSlvr->getInfinity();

     // static part II: B's columns- - - - - - - - - - - - - - - - - - - - - -
     // r[ i ] <= B[ i ] x[ i ] ==> 0 <= B[ i ] x[ i ] - \rho r[ i ] < +INF
     // B[ i ] x[ i ] <= p[ i ] ==> -INF < B[ i ] x[ i ] - \rho p[ i ] <= 0

     if( BNR )  // ... if any
      for( int k = Bbeg[ j ] ; k < Bbeg[ j + 1 ] ; ++k ) {
       // surely one row at position rowpos[ i - 1 ][ Bind[ k ] ]
       Index row = Bind[ k ];
       Index pos = tempI[ count ] = rowpos[ i - 1 ][ row ];

       if( lhs[ row ] && ( lhs[ row ] > -Inf< double >() ) &&
	   rhs[ row ] && ( rhs[ row ] < Inf< double >() ) &&
	   ( lhs[ row ] < rhs[ row ] ) ) {  // two rows are created
        tempHP[ count++ ] = Bval[ k ];  // the first at position pos
        tempI[ count ] = ++pos;         // the second at position pos + 1
        }

       // else only one row is created, but in each case
       tempHP[ count++ ] = Bval[ k ];	// ... with the original value

       }  // end( for( k ) )

     // dynamic part: A's columns- - - - - - - - - - - - - - - - - - - - - - -
     // note that the Lagrangian cost is
     //
     //    c[ i ] - Lambda * A[ i ]
     //
     // which in fact translates into
     //
     //    c[ i ] - ( Lambda + d ) * A[ i ]
     //
     // where Lambda is the stability centre (fixed) and d the direction
     // (the actual variable in the master problem). the fixed term
     // - Lambda * A[ i ] goes into the reduced cose (see below), while the
     // term - d * A[ i ] goes in the coefficient matrix; but when this
     // happens the "-" goes away, see the comments inside ReadFiBLambda()
     // for a complete derivation of the master problem
     for( int k = Abeg[ j ] ; k < Abeg[ j + 1 ] ; ++k ) {
      tempI[ count ] = Aind[ k ] + comp_row[ NrFi ];
      tempHP[ count++ ] = Aval[ k ];
      }

     // add the j-th column of the i-th component- - - - - - - - - - - - - - -
     // note that the objective function of the easy component is
     //
     //    ( c[ i ] - Lambda * A[ i ] ) x[ i ]
     //
     // which in fact translates into
     //
     //    ( c[ i ] - ( Lambda + d ) * A[ i ] ) x[ i ]
     //
     // where Lambda is the stability centre (fixed) and d the direction
     // (the actual variable in the master problem). the cost of the x[ i ]
     // variable then is
     //
     //    c[ i ] - Lambda * A[ i ]
     //
     // clearly, to compute this cost one needs Lambda, but this is not
     // yet available here; hence, THE COST IS CONSTRUCTED ASSUMING THAT
     // Lambda == 0, I.E., AS JUST c[ i ].
     //
     // if the initial stability centre is different from zero, this will
     // have to be explicitly updated by calling ChangeCurrPoint()

     osiSlvr->addCol( count , tempI , tempHP , lbj , ubj , cst[ j ] );

     MSG( 0 , "The variables "<< j << "of the easy component "<< i
              << "has been added" << std::endl );

     }  // end( for each variable j of the component i )

    // now make a final sweep of the rows- - - - - - - - - - - - - - - - - - -
    // this is to put the corresponding (-) RHS into RhoCol

    for( Index j = 0 ; j < BNR ; ++j ) {  // for each row j
     Index pos = rowpos[ i - 1 ][ j ];    // going at position pos

     if( lhs[ j ] && ( lhs[ j ] > -Inf< double >() ) &&
	 rhs[ j ] && ( rhs[ j ] < Inf< double >() ) &&
	 ( lhs[ j ] < rhs[ j ] ) ) {  // two rows are created
      RhoColBse[ RhoColBDm ] = pos;        // first row at position pos
      RhoCol[ RhoColBDm++ ] = - lhs[ j ];  // has nonzero RHS
      RhoColBse[ RhoColBDm ] = ++pos;      // second row at position pos + 1
      RhoCol[ RhoColBDm++ ] = - rhs[ j ];  // has nonzero RHS
      }
     else {                           // one row is created at position pos
      double val = 0;                 // ... but it may have a 0 RHS
      if( lhs[ j ] && ( lhs[ j ] > -Inf< double >() ) )
       val = - lhs[ j ];
      else
       if( rhs[ j ] && ( rhs[ j ] < Inf< double >() ) )
	val = - rhs[ j ];

      if( val ) {
       RhoColBse[ RhoColBDm ] = pos;
       RhoCol[ RhoColBDm++ ] = val;
       }
      }
     }  // end( for each row j )

    // deallocate temporaries- - - - - - - - - - - - - - - - - - - - - - - - -

    delete[] Aval;
    delete[] Aind;
    delete[] Abeg;

    delete[] rhs;
    delete[] lhs;
    delete[] Bval;
    delete[] Bind;
    delete[] Bbeg;

    delete[] ubd;
    delete[] lbd;
    delete[] cst;

    }  // end( if( i is an easy component ) )
   else {  // a difficult component- - - - - - - - - - - - - - - - - - - - - -
           //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // for each component "i" we must create a variable \gamma^i, whose
    // coefficient into objective is - lower bound of the i-th component. the
    // coefficient is 0 if no lower bound is given, as in the rho case.
    // however, unlike in the rho case, \gamma^i must to be fixed to 0
    // whenever some item exists and no bound is provided. \gamma^i lives in
    // the crucial simplex constraint of the i-th component
    //
    //   \sum_k \theta_k^i + \gamma^i = \rho
    //
    // (which is actually a simplex constraint only if \rho is fixed to 1).
    // when the master problem is constructed (now) the bundle is empty and
    // therefore there are no \theta_k^i, which means that \gamma^i is
    // created "free" (>= 0)

    tempI[ 0 ] = comp_row[ i - 1 ];  // the "simplex constraint" for i
    tempHP[ 0 ] = 1;
    osiSlvr->addCol( 1 , tempI , tempHP , 0 , osiSlvr->getInfinity() , 0 );
    MSG( 0 , "The variable gamma_" << i << " has been added and its name is "
         << comp_col[ i ] << std::endl );

    // as outlined above, the "simplex constraint" has a coefficient -1 for
    // \rho: add it now to RhoCol
    RhoColBse[ RhoColBDm ] = comp_row[ i - 1 ];
    RhoCol[ RhoColBDm++ ] = -1;
    }
   }  // end( for( all components ) )

  // after the thusly constructed "static" part there will be the "dynamic"
  // part corresponding to the subgradient, and in particular the vector b
  // of the linear component yb; record where this starts in RhoCol
  RhoColSgPos = RhoColBDm;

  // define the slack and stabilization variables- - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( stab == quadratic ) {
   switchToQP();
   dict_stab = new Index[ MaxSGLen ];
   VectAssign( dict_stab , InINF , MaxSGLen );
   }
  else
   dict_stab = nullptr;

  for( Index i = 0 ; i < CrrSGLen ; ++i ) {
   // ask FiOracle if there are some lower and/or upper bounds on the primal
   bool slack_p = false;       // variables  - - - - - - - - - - - - - - - - -
   bool slack_m = false;

   if( ( Upper[ i ] = Oracle->GetUB( i ) ) < LMINF )
    slack_m = true;  // there is some upper bound on variable i
   if( Oracle->GetUC( i ) )
    Lower[ i ] = -LMINF;
   else {            // the primal variable is constrained to be nonnegative
    Lower[ i ] = 0;
    slack_p = true;
    }

   if( slack_p )
    NNVars++;

   if( slack_p || slack_m )
    BxdVars++;

   // make decision on stabilization - - - - - - - - - - - - - - - - - - - - -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   switch( stab ) {
    case none:
    case boxstep:
     slack_p = slack_m = true;
     break;
    case quadratic:
     tempI[ 0 ] = comp_row[ NrFi ] + i;
     tempHP[ 0 ] = 1;
     dict_stab[ i ] = osiSlvr->getNumCols();
     osiSlvr->addCol( 1 , tempI , tempHP , - osiSlvr->getInfinity() ,
		                             osiSlvr->getInfinity() , 0 );
     MSG( 0 , "The stabilization variable z_" << i
	  << " has been created and its name is " << dict_stab[ i ]
	  << std::endl );
     break;
    default:
     throw( NDOException( "OSIMPSolver::SetDim: undecided stabilization" ) );
    }

   // define the slack variables and put in the vocabulary, one of each pair -
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( slack_p ) {
    tempI[ 0 ] = comp_row[ NrFi ] + i;
    tempHP[ 0 ] = 1;
    osiSlvr->addCol( 1 , tempI , tempHP , 0 , osiSlvr->getInfinity() , 0 );
    dict_slack[ i ] = osiSlvr->getNumCols () - 1;
    MSG( 0 , "The slack s_" << i << "+ has been created and its name is"
    	     << dict_slack[ i ] << std::endl );
    }

   if( slack_m ) {
    tempI[ 0 ] = comp_row[ NrFi ] + i;
    tempHP[ 0 ] = -1;
    osiSlvr->addCol( 1 , tempI , tempHP , 0 , osiSlvr->getInfinity() , 0 );
    if( ! slack_p )
     dict_slack[ i ] = osiSlvr->getNumCols() - 1;
     MSG( 0 , "The slack s_" << i << "- has been created and its name is "<<
    	( slack_p? dict_slack[ i ] + 1 : dict_slack[ i ] ) << std::endl );
     }

   }  // end( for( all variable i ) )

  // add the column of rho in the master - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  comp_col[ 0 ] = osiSlvr->getNumCols();
  osiSlvr->addCol( RhoColBDm , RhoColBse , RhoCol , 0 ,
		   osiSlvr->getInfinity() , 0 );
  MSG( 0 , "rho has been created and its name is " << comp_col[ 0 ]
           << std::endl );

  // change the price of the slack variables - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tUpdatePrices();       // change the prices for the stabilization
  if( stab != boxstep )  // also change the slack prices (if they have not
   ptUpdatePrices();     // already been changed in tUpdatePrices())

  // final memory allocations- - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // allocate the memory for the lower bounds
  LwrBnds = new HpNum[ NrFi + 1 ];
  VectAssign( LwrBnds , -HpINF , NrFi + 1 );  // Lower Bounds == none

  // allocate the memory for Gi * d
  GiPerd = new HpNum[ MaxBSize ];

  }  // end( if( Oracle ) )
 else {  // sets the max bundle size to n and activate/deactivate the Active
         // Set Mechanism without changing anything else - - - - - - - - - - -

  if( MxBSz != MaxBSize ) {
   // delete the item in excess- - - - - - - - - - - - - - - - - - - - - - - -
   // the existing items in the bundle (if any) are all kept if MxBSz is
   // larger than MaxBSize, but a smaller value will force deletion of all
   // the items with "name"  >= MxBSz

   if( MxBSz < MaxBSize ) {
    // eliminate all items with name >= new size (keep only names < new size)
    // remember that item_maxname is "maximum name + 1"
    Index mxnm = item_maxname;
    for( ; --mxnm >= MxBSz ; )
     if( dict_item[ mxnm ] != InINF )
      RmvItem( mxnm );

    // find if the maximum name is smaller yet
    while( dict_item[ mxnm ] == InINF )
     mxnm--;

    item_maxname = mxnm + 1;
    }

   // resize the data structures - - - - - - - - - - - - - - - - - - - - - - -
   Index_Set olddict_item = dict_item;
   Index_Set oldwcomp = wcomp;
   HpRow oldGiPerd = GiPerd;
   
   dict_item = new Index[ MxBSz ];
   wcomp = new Index[ MxBSz ];
   GiPerd = new HpNum[ MxBSz ];

   VectAssign( dict_item , InINF , MxBSz );
   VectAssign( wcomp , InINF , MxBSz );

   cIndex copydim = ( MxBSz < MaxBSize ) ? MxBSz : MaxBSize;

   if( olddict_item ) {
    VectAssign( dict_item , olddict_item , copydim );
    delete[] olddict_item;
    }

   if( oldwcomp ) {
    VectAssign( wcomp , oldwcomp , copydim );
    delete[] oldwcomp;
    }

   if( oldGiPerd ) {
    VectAssign( GiPerd , oldGiPerd , copydim );
    delete[] oldGiPerd;
    }
   
   MaxBSize = MxBSz;  // new dimension of the bundle

   // the initial "active set" of variables is *empty* and so all the
   // variables must be non active - - - - - - - - - - - - - - - - - - - - - -
   if( ( ! useactiveset ) && UsAvSt )
    for( Index i = 0 ; i < CrrSGLen ; i++ )
     deactivate( i );

   // all the  variables are considered to be always "active"- - - - - - - - -
   if( useactiveset && ( ! UsAvSt ) )
    for( Index i = 0 ; i < CrrSGLen ; i++ )
     activate( i );

   useactiveset = UsAvSt;

   }  // end( if( MxBSz != MaxBSize ) )
  }  // end( else( Oracle ) )- - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 }  // end( OSIMPSolver::SetDim )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::Sett( cHpNum tt )
{
 if( t != tt ) {  // actually update t- - - - - - - - - - - - - - - - - - - -
  t = tt;
  MSG( 0 , "OSIMPSolver::Sett(): t = " << t << std::endl );
  tUpdatePrices();
  }
 }  // end( OSIMPSolver::Sett )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetPar( const int wp , cHpNum value )
{
 switch( wp ) {
  case( kMaxTme ): MaxTime = value; break;

  case( kOptEps ):  // because the dual formulation is used, the role of
   OptEps = value;  // the primal and dual tolerances is reversed
   osiSlvr->setDblParam( OsiPrimalTolerance , double( value ) );
   break;

  case( kFsbEps ):  // because the dual formulation is used, the role of
   FsbEps = value;  // the primal and dual tolerances is reversed
   osiSlvr->setDblParam( OsiDualTolerance , double( FsbEps ) );
   #if WHICH_OSI_MP == 1
    if( algorithm == kBar ) {
     auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
     CPXsetdblparam( osiCpx->getEnvironmentPtr() ,
		     CPX_PARAM_BAREPCOMP , value );
     }
   #endif
   break;

   // case( kZero ): break;  // just ignore it- - - - - - - - - - - - - - - -

  default: throw( NDOException( "SetPar( HpNum ): unknown parameter" ) );
  }
 }  // end( OSIMPSolver::SetPar( HpNum ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetThreads( int nthreads )
{
 if( ! osiSlvr )
  return;

 #if WHICH_OSI_MP == 1
  auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
  CPXsetintparam( osiCpx->getEnvironmentPtr() ,
		  CPXPARAM_Threads , nthreads );
 #elif WHICH_OSI_MP == 2
  auto osiGrb = dynamic_cast< OsiGrbSolverInterface * >( osiSlvr );
  GRBsetintparam( osiGrb->getEnvironmentPtr() , "Threads" , nthreads );
 #endif

 
 }  // end( OSIMPSolver::SetThreads )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetLowerBound( cHpNum LwBnd , cIndex wFi )
{
 // important note: each "hard" component has a "slack" variable \gamma_i in
 // the corresponding convexity constraint, and there is a single "slack"
 // variable \rho that enters *all* the convexity constraints. \rho is used to
 // represent the global lower bound (its cost is - that, where the - comes
 // from the fact that r + \rho = 1, where r would be the slack variable):
 // when there is no global lower bound, \rho is fixed to 0 (by fixing its
 // upper bound). the \gamma_i are similar in that their costs are the
 // individual lower bounds (without the "-"), i.e., the *opposite* of the
 // (usually, non-negative) linearization error of the all-0 "flat"
 // subgradient corresponding to the lower bound. if there is no individual
 // lower bound for component i, \gamma_i should be fixed to 0 (hence its
 // cost is irrelevant). however, \gamma_i also serves another purpose: if
 // there are no subgradients of component i, which would make it impossible
 // to satisfy the convexity constraint making the master problem empty
 // (note that constraints do no impact the convexity constraint anyway and
 // therefore cannot help), \gamma_i is made free even if there is no lower
 // bound to make the master problem feasible. when this is done, it is not
 // exactly trivial to chose its cost, which corresponds to choosing the a
 // "fake" lower bound for the model of component i.
 //
 // the primal of the solved model contains the constraints
 //
 //     v >= \sum_{ all components } < value of model for that component >
 //
 //     v >= l
 //
 // where v is the total (translated) model value and l is the (translated)
 // global lower bound (actually something to this effect due to the
 // reformulation of the r + \rho = 1 constraint, see comment in
 // ReadFiBLambda() for details). note that only "hard" components have
 // lower bounds ("easy" ones have none for obvious reasons) and the cost of
 // \gamma_i for "hard" component i is the lower bound l_i on the
 // corresponding variable v_i representing the value of the (translated)
 // model in the above constraint. clearly, when the bound l_i is "fake",
 // i.e, the variable is free only to make the problem feasible but there is
 // no actual finite lower bound on component i, the value of l_i has to be
 // chosen in such a way that the first constraint do not mistakenly become
 // more stringent than the second. this is done by chosing l_i < l.
 
 MSG( 0 , "OSIMPSolver::SetLowerBound()\n" );

 if( ( wFi < InINF ) && weasy[ wFi ] )
  throw( NDOException( "OSIMPSolver::SetLowerBound: wFi easy component" ) );

 cIndex h = wFi > NrFi ? 0 : wFi;  
 if( LwrBnds[ h ] == LwBnd )  // setting an already known bound
  return;                     // nothing to do

 cIndex cc = comp_col[ h ];   // the column corresponding to \gamma_i

 if( LwBnd > -HpINF ) {  // setting the lower bound- - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( h ) {  // an idividual lower bound
   // if the lower bound was not set already, "free" \gamma_i by setting its
   // upper bound to +INF; however this has to be done only if there are
   // subgradients of component i in the bundle, for otherwise \gamma_i is
   // free already
   if( NSubG[ h ] && ( LwrBnds[ h ] == -HpINF ) )
    osiSlvr->setColUpper( cc , osiSlvr->getInfinity() );

   // set LwBnd as coefficient of \gamma_i in the objective function
   osiSlvr->setObjCoeff( cc , LwBnd );

   MSG( 0 , "LB( " << wFi << " ) = " << LwBnd << std::endl );
   }
  else {     // the global lower bound
   // note that the constraint on \rho needs change:
   // r + \rho = 1 ----> r = 1 - \rho => 0 [ \rho <= 1 ]
   // and in the objective function LwBnd * r = LwBnd - LwBnd * \rho
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   osiSlvr->setRowBounds( 0 , - osiSlvr->getInfinity() , 1 );
   osiSlvr->setObjCoeff( cc , - LwBnd );

   // for those non-easy components that currently have neither subgradients
   // nor individual bound, ensure that the cost of the corresponding
   // \gamma_i is < than (-) that of \rho
   for( Index i = 0 ; ++i <= NrFi ; )
    if( ( ! weasy[ i ] ) && ( ! NSubG[ i ] ) && ( LwrBnds[ i ] == -HpINF ) )
     osiSlvr->setObjCoeff( comp_col[ i ] , LwBnd - 1 );
     
   MSG( 0 , "Total LB = " << LwBnd << std::endl );
   }
  }
 else {  // LwBnd == -HpINF: reset the lower bound - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( h ) {    // an individual lower bound

   // fix \gamma_h to 0, unless there are no subgradients of component h in
   // the bundle, in which case it has to be kept free
   if( NSubG[ h ] ) {
    osiSlvr->setColUpper( cc , 0 );
    osiSlvr->setObjCoeff( cc , 0 );
    }
   else {
    // if \gamma_h is left free, which means there are no subgradients of
    // component h in the bundle, its cost depends on whether or not a
    // global lower bound is set: if not it is 0 and thus remains, but if
    // so it has to be set to < than (-) that of \rho
    if( LwrBnds[ 0 ] > -HpINF )
     osiSlvr->setObjCoeff( cc , LwrBnds[ 0 ] - 1 );
    }

   MSG( 0 , "Lower bound of "<< wFi << " has been reset " << std::endl );
   }
  else {       // the global lower bound
   // change the \rho constraint --> \rho = 1 and put the coefficient
   // of \rho equal to zero in the objective function

   MSG( 0 , "Total lower bound is - Inf< double >()\n" );
   osiSlvr->setObjCoeff( cc , 0 );
   osiSlvr->setRowBounds( 0 ,  1 , 1 );

   // for those non-easy components that currently have neither subgradients
   // nor individual bound, and whose objective was set to a "strange"
   // value to avoid interfacing with \rho, just put a 0 there
   for( Index i = 0 ; ++i <= NrFi ; )
    if( ( ! weasy[ i ] ) && ( ! NSubG[ i ] ) &&	( LwrBnds[ i ] == -HpINF ) )
     osiSlvr->setObjCoeff( comp_col[ i ] , 0 );
 
   MSG( 0 , "Global lower bound has been reset " << std::endl );
   }
  }

 LwrBnds[ h ] = LwBnd;  // finally, record the lower bound

 }  // end( OSIMPSolver::SetLowerBound )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetMPLog( std::ostream * outs , const char lvl )
{
 MPSolver::SetMPLog( outs , lvl );

 #if( WHICH_OSI_MP == 1 )
  auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
  auto env = osiCpx->getEnvironmentPtr();

   #if( OSIMPSOLVERLOG )
    if( outs && MPLLvl ) {
     CPXsetlogfilename( env , "cplex.log" , "w" );
     CPXsetintparam( env , CPXPARAM_MIP_Display , MPLLvl );
     }
    else
   #endif
   {
    CPXsetlogfilename( env , NULL , NULL ) ;

    CPXsetintparam( env , CPX_PARAM_SCRIND , CPX_OFF );
    CPXsetintparam( env , CPXPARAM_MIP_Display , CPX_OFF );
    CPXsetintparam( env , CPXPARAM_Advance , 0 );
    CPXsetintparam( env , CPXPARAM_ScreenOutput , CPX_OFF );
    CPXsetintparam( env , CPXPARAM_Tune_Display , CPX_OFF );
    CPXsetintparam( env , CPXPARAM_Barrier_Display , 0 );
    CPXsetintparam( env , CPXPARAM_Simplex_Display , 0 );
    CPXsetintparam( env , CPXPARAM_Sifting_Display , 0 );
    CPXsetintparam( env , CPXPARAM_Network_Display , 0 );
    CPXsetintparam( env , CPXPARAM_ParamDisplay , CPX_OFF );
    }
 #elif( WHICH_OSI_MP == 2 )
  auto osiGrb = dynamic_cast< OsiGrbSolverInterface * >( osiSlvr );
  auto env = osiGrb->getEnvironmentPtr();

  #if( OSIMPSOLVERLOG )
   if( outs && MPLLvl )
    GRBsetstrparam( env , "LogFile" , "gurobi.log" );
   else
    GRBsetstrparam( env , "LogFile" , "" );

   GRBsetintparam( env , "OutputFlag" , MPLLvl );
  #else
   GRBsetintparam( env , "OutputFlag" , 0 );
  #endif
 #endif

 } // end( OSIMPSolver::SetMPLog )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetAlgo( const OsiAlg algo , const OsiRed reduc )
{
 algorithm = algo;

 #if WHICH_OSI_MP == 1
  // setting algorithmic parameters for Cplex- - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
  auto env = osiCpx->getEnvironmentPtr();

  // primal and dual reduction - - - - - - - - - - - - - - - - - - - - - - -
  switch ( reduc ) {
   case( rNo ):
    CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_NOPRIMALORDUAL );
    break;
   case( rPrim ):
    CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_PRIMALONLY );
    break;
   case( rDual ):
    CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_DUALONLY );
    break;
   case( rBoth ):
    CPXsetintparam( env , CPX_PARAM_REDUCE , CPX_PREREDUCE_PRIMALANDDUAL );
   }

 if( stab != quadratic )  // choose CPLEX algorithm- - - - - - - - - - - - -
  switch( algorithm ) {
   case( kAuto ): CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_AUTOMATIC );
                  break;
   case( kPrim ): CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_PRIMAL );
                  break;
   case( kDual ): CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_DUAL );
                  break;
   case( kNet ):  CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_NET );
                  break;
   case( kBar ):  CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_BARRIER );
                  break;
   case( kSif ):  CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_SIFTING );
                  break;
   case( kCon ):  CPXsetintparam( env , CPX_PARAM_LPMETHOD ,
				  CPX_ALG_CONCURRENT );
   } // end switch
 else
  switch( algorithm ) {
   case( kAuto ): CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_AUTOMATIC );
                  break;
   case( kPrim ): CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_PRIMAL );
                  break;
   case( kDual ): CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_DUAL );
                  break;
   case( kNet ):  CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_NET );
                  break;
   case( kBar ):  CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_BARRIER );
                  break;
   case( kSif ):  CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_SIFTING );
                  break;
   case( kCon ):  CPXsetintparam( env , CPX_PARAM_QPMETHOD ,
				  CPX_ALG_CONCURRENT );
   } // end switch
 #elif WHICH_OSI_MP == 2
  // setting algorithmic parameters for Gurobi - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  auto osiGrb = dynamic_cast< OsiGrbSolverInterface * >( osiSlvr );
  auto env = osiGrb->getEnvironmentPtr();

  // primal and dual reduction - - - - - - - - - - - - - - - - - - - - - - -
  if( ( reduc == rNo ) || ( reduc == rDual ) )
   GRBsetintparam( env , "Presolve" , int( 0 ) );
  else
   GRBsetintparam( env , "Presolve" , int( -1 ) );
   
  if( ( reduc == rNo ) || ( reduc == rPrim ) )
   GRBsetintparam( env , "DualReductions" , int( 0 ) );
  else
   GRBsetintparam( env , "DualReductions" , int( 1 ) );

  switch( algorithm ) {
   case( kAuto ): GRBsetintparam( env , "Method" , int( -1 ) );
                  GRBsetintparam( env , "Sifting" , int( -1 ) );
                  break;
   case( kPrim ): GRBsetintparam( env , "Method" , int( 0 ) );
                  GRBsetintparam( env , "Sifting" , int( -1 ) );
                  break;
   case( kDual ): GRBsetintparam( env , "Method" , int( 1 ) );
                  GRBsetintparam( env , "Sifting" , int( -1 ) );
                  break;
   case( kNet ):  throw( NDOException(
			   "OSIMPSolver::SetAlgo: Gurobi has no Net" ) );
                  break;
   case( kBar ):  GRBsetintparam( env , "Method" , int( 2 ) );
                  GRBsetintparam( env , "Sifting" , int( -1 ) );
                  break;
   case( kSif ):  GRBsetintparam( env , "Method" , int( -1 ) );
                  GRBsetintparam( env , "Sifting" , int( 2 ) );
                  break;
   case( kCon ):  GRBsetintparam( env , "Method" , int( 3 ) );
                  GRBsetintparam( env , "Sifting" , int( -1 ) );
   } // end switch
 #endif

 } // end ( OSIMPSolver::SetAlgo )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetOsi( OsiSolverInterface *osi )
{
 if( osiSlvr && osi && FIO )
  throw( NDOException(
              "OSIMPSolver::SetOsi: cannot change OSI solver in flight" ) );
 delete osiSlvr;
 osiSlvr = osi;

 if( osiSlvr ) {
  #if WHICH_OSI_MP == 1
   if( ! dynamic_cast< OsiCpxSolverInterface * >( osiSlvr ) )
    throw( NDOException( "OSIMPSolver::SetOsi: not an OsiCpx" ) );
  #elif WHICH_OSI_MP == 2
   if( ! dynamic_cast< OsiGrbSolverInterface * >( osiSlvr ) )
    throw( NDOException( "OSIMPSolver::SetOsi: not an OsiGrb" ) );
  #endif

  // set to 0 the log: do that to avoid some warnings in cleanup
  // osiSlvr->messageHandler()->setLogLevel( int( MPLLvl ) );
  osiSlvr->passInMessageHandler( derhand );

  // initialise the log, likely to MPLLvl == 0; this is useful to avoid
  // that certain underlying solvers start spewing out unwanted stuff
  OSIMPSolver::SetMPLog( 0 , MPLLvl );

  // ensure the solver knows of the tolerances
  SetPar( kOptEps , OptEps );
  SetPar( kFsbEps , FsbEps );
  }
 else
  cleanup();

 }  // end( OSIMPSolver::SetOsi )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetStabType( const StabFun sf )
{
 MSG( 0 , "OSIMPSolver::SetStabType(()\n" );
 if( ( ( stab = sf ) != unset ) && FIO )
  throw( NDOException( 
                 "OSIMPSolver::SetStabType: cannot change D_t in flight" ) );

 }  // end( OSIMPSolver::SetStabType )

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

OSIMPSolver::MPStatus OSIMPSolver::SolveMP( void )
{
 MSG( 0 , "OSIMPSolver::SolveMP()\n" );

 MPStatus MPsts; // status

 if( MPt )
  MPt->Start();

 ++NRCall;

 // get the dimension of the Master Problem  - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 int nc = osiSlvr->getNumCols();
 int nr = osiSlvr->getNumRows();

 #if( OSIMPSOLVERLOG )
  if( osiSlvr && ( MPLLvl > 1 ) && nc && nr ) {
   std::string fname( "MP_" + std::to_string( NRCall ) + ".lp" );
   // use native solver API when possible, because OSI ignores (and,
   // therefore, does not print) the quadratic objective, if any
   #if WHICH_OSI_MP == 1
    auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
    CPXwriteprob( osiCpx->getEnvironmentPtr() , osiCpx->getLpPtr() ,
		  fname.c_str() , "lp" );
   #elif WHICH_OSI_MP == 2
    auto osiGrb = dynamic_cast< OsiGrbSolverInterface * >( osiSlvr );
    GRBwrite( osiGrb->getLpPtr() , fname.c_str() );
   #else
    osiSlvr->writeLp( fname.c_str() );
   #endif
   }
 #endif

 // set algorithmic parameters - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MaxTime ) {
  // amazingly, there is no method in OsiSolverInterface to do this!
  #if WHICH_OSI_MP == 1
   auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
   CPXsetdblparam( osiCpx->getEnvironmentPtr() ,
		   CPX_PARAM_TILIM , MaxTime );
  #elif WHICH_OSI_MP == 2
   auto osiGrb = dynamic_cast< OsiGrbSolverInterface * >( osiSlvr );
   GRBsetdblparam( osiGrb->getEnvironmentPtr() ,
		   GRB_DBL_PAR_TIMELIMIT , MaxTime );
  #else
   throw( NDOException( "OSIMPSolver::SolveMP: cannot set time limit" ) );
  #endif
  }
 
 // solve the Master Problem - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( stab == quadratic ) {
  #if WHICH_OSI_MP == 1
   // call Cplex QP solver - - - - - - - - - - - - - - - - - - - - - - - - -
   auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
   CPXqpopt( osiCpx->getEnvironmentPtr() , osiCpx->getLpPtr() );
  #elif WHICH_OSI_MP == 2
   // call Gurobi QP solver- - - - - - - - - - - - - - - - - - - - - - - - -
   auto osiGrb = dynamic_cast< OsiGrbSolverInterface * >( osiSlvr );
   GRBoptimize( osiGrb->getLpPtr() );
  #else
   throw( NDOException( "OSIMPSolver::SolveMP: not implemented yet" ) );
  #endif
  }
 else
  if( first ) {
   osiSlvr->initialSolve();
   first = false;
   }
  else
   osiSlvr->resolve();

 // check the status of the solver and take the necessary action - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 MPsts = kOK;

 #if( PRESERVE_OSI_SOLS )
  // note: do this independently from the return code of the MPSolver to
  //       ensure that csol and rcst are well-defined even under error

  if( nc > csols ) {
   delete[] csol;
   delete[] rcst;
   csols = nc;
   csol = new double[ csols ];
   rcst = new double[ csols ];

   // csol[ comp_col[ 0 ] ] is checked by ReadMult(), so ensure it it
   // properly initialized in any case
   csol[ comp_col[ 0 ] ] = 0;
   }
 #endif

 if( osiSlvr->isProvenDualInfeasible() ) {  // the feasible set is empty
  MSG( 0 , "MP is dual infeasible" << std::endl );
  MPsts = kUnfsbl;
  }
 else
  if( osiSlvr->isProvenPrimalInfeasible() ) {
   MSG( 0 , "MP is primal infeasible" << std::endl );
   MPsts = kUnbndd;
   }
  else
   if( ! osiSlvr->isProvenOptimal() ) {
    #if WHICH_OSI_MP == 1
     auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
     auto env = osiCpx->getEnvironmentPtr();
     auto qp = osiCpx->getLpPtr();
     if( CPXgetstat( env , qp ) == CPX_STAT_ABORT_TIME_LIM ) {
      MSG( 0 , "MP stopped by time limit" << std::endl );
      MPsts = kStppd;
      }
    #elif WHICH_OSI_MP == 2
     auto osiGrb = dynamic_cast< OsiGrbSolverInterface * >( osiSlvr );
     int status;
     GRBgetintattr( osiGrb->getLpPtr() , "Status" , &status );
     if( status == 9 ) {
      MSG( 0 , "MP stopped by time limit" << std::endl );
      MPsts = kStppd;
      }
     #endif

    if( osiSlvr->isAbandoned() ) {
     MSG( 0 , "Warning: numerical difficulties in the solver, ignoring them"
	   << std::endl );
     }
    else {
     MSG( 0 , "Some unknown error happened" << std::endl );
     MPsts = kError;
     }
    } // end status check - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // record solution information in its data structures - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MPsts == kOK ) {

  #if( PRESERVE_OSI_SOLS )
   // note: for the quadratic MP one can still call getColSolution(),
   // getReducedCost() and getRowPrice(), as Cplex uses the same functions
   // as in the LP case (and precisely CPXgetx, CPXgetdj and CPXgetpi)

   VectAssign( csol , osiSlvr->getColSolution() , nc );
   VectAssign( rcst , osiSlvr->getReducedCost() , nc );

   if( nr > rsols ) {
    delete[] rsol;
    rsols = nr;
    rsol = new double[ rsols ];
    }

   VectAssign( rsol , osiSlvr->getRowPrice() , nr );
  #else
   const double *rsol = osiSlvr->getRowPrice();
   const double *rcst = osiSlvr->getReducedCost();
  #endif

  // to compute Gid use the reduced cost rc_i, in particular we have:
  // g^{top} * delta = rc_i + v_i + \alpha_i

  const double *objcoeff = osiSlvr->getObjCoefficients();
  for( Index i = 0 ; i < item_maxname ; i++ )
   if( dict_item[ i ] < InINF ) {
    GiPerd[ i ] = rcst[ dict_item[ i ] ] - objcoeff[ dict_item[ i ] ];
    if( IsSubG( i ) )
     GiPerd[ i ] += rsol[ comp_row[ WComp( i ) - 1 ] ];
    }

  // print the primal/dual solution  - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( OSIMPSOLVERLOG )
   if( MPLLvl > 2 ) {
    #if( PRESERVE_OSI_SOLS == 0 )
     const double *csol = osiSlvr->getColSolution();
    #endif

    *MPLog << "x =" << std::endl;
    for( int i = 0 ; i < nc ; i++ )
     *MPLog << csol[ i ]  << " ~ ";
    *MPLog << std::endl;

    *MPLog << "Pi =" << std::endl;
    for( int i = 0 ; i < nr ; i++ )
     *MPLog << rsol[ i ] << " ~ ";
    *MPLog << std::endl;
    }
  #endif

  // mark all items as read  - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < item_maxname ; i++ )
   if( wcomp[ i ] < InINF )
    wcomp[ i ] &= ~LLBIndex;

  } // end if - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( MPt )
  MPt->Stop();

 return( MPsts );

 }  // end( OSIMPSolver::SolveMP )

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR READING RESULTS ------------------------*/
/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadFiBLambda( cIndex wFi )
{
 // for justifying some of the formulae below, let us start by recalling
 // that for a "partial quadratic" problem
 //
 //  min { q' x' + q" x" + x"T Q x" / 2 : A' x' + A" x" >= b }
 //
 // the Lagrangian is
 //
 //  min { q' x' + q" x" + x"T Q x" / 2 + y ( b - A' x' + A" x" ) }
 //
 // and thus the Lagrangian Dual is
 //
 //  max { < the above > : y >= 0 }
 //
 // rearranging this gives
 //
 //  max { y b + min{ ( q' - y A' ) x' : x' }
 //            + min { ( q" - y A" ) x" + x"T Q x" / 2 : x" } }
 //
 // the first subproblem is only finite when y A' = q'
 //
 // the second subproblem's optimality condition is
 //
 //  ( q" - y A" ) + Q x" = 0  ==>  x" = Q^{-1} ( y A" - q" )
 //
 // substituting in the objective gives
 //
 //  ( q" - y A" ) Q^{-1} ( y A" - q" ) +
 //  ( y A" - q" ) Q^{-1} Q Q^{-1} ( y A" - q" ) / 2 
 //  = - ( y A" - q" ) Q^{-1} ( y A" - q" ) / 2 
 //
 // this is better rewritten with auxiliary variables
 //
 //  - s^T Q^{-1} s / 2 : s = y A" - q" 
 //
 // finally yielding the "partial quadratic dual"
 //
 //  max { y b - s^T Q^{-1} s / 2 : y A' = q' , s = y A" - q" , y >= 0 }
 
 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif

 HpNum FiVal = 0;
 if( ! wFi ) {  // the 0-th component- - - - - - - - - - - - - - - - - - - -
  // it is given by Fi[ 0 ]( d^* ) = d^* b. Note that d^* is the dual
  // solution of the dynamic part and it represents the search direction

  for( Index i = RhoColSgPos ; i < RhoColBDm ; ++i )
   FiVal -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];

  MSG( 0 , "FiBLambda( 0 ) = " << FiVal << std::endl );
  }
 else
  if( wFi <= NrFi ) {  // a specific component - - - - - - - - - - - - - - -
   // the easy case is quite different from the difficult one: in fact, in
   // the easy case it computes the actual value of the function while in
   // the other case just an approximation is provided

   if( weasy[ wFi ] ) {
    // compute the actual function value of the easy components by using the
    // row prices. here is the full development of the formula
    //
    // we start with the simple version of the problem
    //
    //  min { f(x) + max { ( c - xA ) u : B u <= b } }
    //
    // with one hard component an an easy one written in the most compact
    // possible form (hiding the complexity due to two-sided constraints
    // and variable bounds that possibly result in duplicated rows). the
    // corresponding stabilised primal master problem reads
    //
    //  min  v + max { ( c - ( x + d ) A ) u : B u <= b } + 1/(2t) d^T d 
    //       v e >= G d - a
    //
    // where we hide the dependency of a from the stability centre x. by
    // dualising the inner max we get
    //
    //  min { w b : w B = c - ( x + d ) A , w >= 0 }
    //
    // and hence
    //
    //  min    v +  w b + 1/(2t) d^T d 
    //       e v             - G d      >= - a           (y)
    //              w B      +   d A    == c - x A = c'  (u)
    //              w >= 0
    //
    // whose dual (which is what we ideally solve) therefore is
    //
    //  max - y a + c' u - t s^T s / 2
    //        y e                         == 1           (v)
    //               B u                  <= b           (w)
    //      - y G  + A u   + s            == 0
    //        y >= 0
    //
    // from the complementary slackness conditions
    //
    //     w ( b - B u ) = 0  ==> w b = w B u
    //
    // and from
    //
    //      ( w B + d A - c + x A ) u = 0
    //
    // we get
    //
    //      w b = w B u = ( c - ( x + d ) A ) u
    //
    // hence, the value of the easy component in the optimal primal solution
    // d (corresponding to the optimal dual solution u) can be computed in
    // linear time as just w b
    //
    // however, in practice things are more complex because of the global
    // lower bound. that is, the primal problem actually is
    //
    //  min  v              + 1/(2t) d^T d 
    //       v -   v' - w b                >= 0             (z)
    //       v                             >= l             (r)
    //           e v'            - G d     >= - a           (y)
    //                  w B      +   d A   == c - x A = c'  (u)
    //                  w >= 0
    //
    // and hence the dual is
    //
    //  max         r l - y a + c' u - t s^T s / 2
    //        z   + r                                == 1   (v)
    //      - z         + y e                        == 0   (v')
    //      - z b              + B u                 <= 0   (w)
    //                  - y G  + A u   + s           == 0
    //        z >= 0
    //              r >= 0
    //                    y >= 0
    //
    // thus, the relevant complementary slackness now are
    //
    //     w ( z b - B u ) = 0  ==> z ( w b ) = w B u
    //
    // which together with the unchanged
    //
    //      ( w B + d A - c + x A ) u = 0
    //
    // gives
    //
    //      z ( w b ) = w B u = ( c - ( x + d ) A ) u
    //
    // i.e., finally the value of the easy component as
    //
    //       ( c - ( x + d ) A ) u = z ( w b )
    //
    // where z == comp_col[ 0 ]. thus, if a global lower bound is defined,
    // w b has to be scaled by z IF ONE WANTS THE CONTRIBUTION OF THE EASY
    // COMPONENT TO THE VALUE OF v. however, THIS IS NOT WHAT WE RETURN
    // HERE, WHICH RATHER IS THE "TRUE" VALUE OF THE "EASY" COMPONENT
    // IRRESPEVTIVE OF THE LOWER BOUND, which is just the un-scaled w b.
    //
    // note that the same distinction would hold for "hard" components: the
    // contribution of the "hard" component to v is v' only if there is no
    // global lower bound l, otherwise it should be scled by z as well. but,
    // as for "easy" components, this is not the value that is returned here,
    // since the un-scaled one is rather the useful one.
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // all the entries of the RHS corresponding to the component wFi are
    // consecutive, so seek for the first one exploiting the fact that
    // RhoColBse is ordered
    const auto stp = RhoColBse + RhoColSgPos;
    auto it = std::lower_bound( RhoColBse , stp , comp_row[ wFi - 1 ] );
    for( auto rit = RhoCol + std::distance( RhoColBse , it ) ; it < stp ; ) {
     Index pos = *(it++);
     if( pos >= comp_row[ wFi ] )
      break;
     FiVal -= rsol[ pos ] * (*(rit++));
     }

    /* if one wanted the contribution of the easy component to the value of
       v, which we don't, it would be necessary to scale by z
    if( LwrBnds[ 0 ] > -HpINF )
     FiVal *= csol[ comp_col[ 0 ] ];
    */
    }
   else
    // for "hard" components, the value is significant only if there is at
    // least a subgradient in the bundle, comprised the "implicit all-0 one"
    // corresponding to a finite lower bound
    if( NSubG[ wFi ] || ( LwrBnds[ wFi ] >- HpINF ) )
     FiVal = rsol[ comp_row[ wFi - 1 ] ];
    else
     FiVal = HpINF;

   MSG( 0 , "FiBLambda(" << wFi << ") = " << FiVal << std::endl );
   }
  else { // all the components - - - - - - - - - - - - - - - - - - - - - - -
   // the value of v*, the global model, is the dual of the r + \rho = 1
   // constraint. however, when the lower bound is added that constraint
   // is transformed into \rho <= 1, and therefore the dual need to be
   // changed. here is the full development of the formula
   //
   // we start with a stylised version of the master problem with two hard
   // components, no easy one, and the individual lower bounds "hidden"
   // among the original subgradients:
   //
   //  min  v               + 1/(2t) d^T d 
   //       v    - v1    - v2               >= 0             (z)
   //       v                               >= l             (r)
   //           e1 v1           - G1 d      >= - a1          (y1)
   //                   e2 v2   - G2 d      >= - a2          (y2)
   //
   // its dual therefore is
   //
   //  max      r l - y1 a1 - y2 a2 - t s^T s / 2
   //       z + r                                 = 1        (v)
   //     - z       + y1 e1                       = 0        (v1)
   //     - z               + y2 e2               = 0        (v2)
   //               - y1 G1 - y2 G2        - s    = 0        (d)
   //       z >= 0 , r >= 0 , y1 >= 0 , y2 >= 0
   //
   // however, this is not exactly the problem we solve, because we do the
   // substitution
   //
   //  r = 1 - z
   //
   // which, considering the sign constraint r >= 0 ==> z <= 1, yields
   //
   //  l +
   //  max - z l - y1 a1 - y2 a2 - t s^T s / 2
   //        z                                  <= 1        (v')
   //      - z   + y1 e1                         = 0        (v1)
   //      - z           + y2 e2                 = 0        (v2)
   //            - y1 G1 - y2 G2         - s     = 0        (d)
   //        z >= 0 , y1 >= 0 , y2 >= 0
   //
   // working backward from that, the dual is
   //
   //  l +
   //  min  v'               + 1/(2t) d^T d 
   //       v'    - v1   - v2               >= - l          (z)
   //       v'                              >= 0            [none, sign]
   //            e1 v1          - G1 d      >= - a1         (y1)
   //                   e2 v2   - G2 d      >= - a2         (y2)
   //
   // now, it is immediate to realise that this problem is completely
   // equivalent to the initial one under the substitution
   //
   //  v = v' + l    <==>     v' = v - l
   //
   // which, starting from
   //
   //  v' + l >= v1 + v2   ,   v' >= 0
   //
   // gives back the original
   //
   //  v >= v1 + v2        ,   v >= l
   //
   // finally justifying the formula below (l is the opposite of the current
   // cost of z == \rho, i.e., comp_col[ 0 ])

   FiVal = rsol[ 0 ];

   if( LwrBnds[ 0 ] > - HpINF )
    FiVal += LwrBnds[ 0 ];
   else {
    // if there is no global lower bound, the problem "implicitly" still
    // has the original constraint
    //
    //       z + r                                 = 1        (v)
    //
    // except that it also implicitly has "r = 0", hence v is directly the
    // dual variable of that constraint (the very first one)

    // if there is no global bound, and some of the "hard" components is
    // undefined, then the total value of the model also is
    for( Index i = 0 ; ++i <= NrFi ; )
     if( ( ! weasy[ i ] ) && ( ! NSubG[ i ] ) &&
	 ( LwrBnds[ i ] == -HpINF ) ) {
      FiVal = HpINF;
      break;
      }
    }

   if( wFi < InINF ) {     // the 0-th component has to be excluded
    if( FiVal < HpINF ) {  // excluded from undefined is undefined
     for( Index i = RhoColSgPos ; i < RhoColBDm ; ++i )
      FiVal += rsol[ RhoColBse[ i ] ] * RhoCol[ i ];
     }

    MSG( 0 , "FiBLambda( all - 0 ) = " << FiVal << std::endl );
    }
   else
    MSG( 0 , "FiBLambda = " << FiVal << std::endl );

   }  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 return( FiVal );

 }  // end( OSIMPSolver::ReadFiBLambda )

/*------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadDt( cHpNum tt )
{
 HpNum value = 0;
 HpNum mytt;

 #if( PRESERVE_OSI_SOLS == 1 )
  cLMRow primal = rsol + comp_row[ NrFi ];
 #else
  cLMRow primal = osiSlvr->getRowPrice() + comp_row[ NrFi ];
  const double *csol = osiSlvr->getColSolution();
 #endif

 switch( stab ) {
  case none: break;
  case boxstep:
  mytt = tt * ( 1 + FsbEps );  // enlarge tt a little to take in account
                               // for numerical errors (use the tolerance)
  if( mytt < t ) {  // check whether d lies in the ball of radious tt
   for( Index n = 0 ; n < CrrSGLen ; n++ )
    if( ABS( primal[ n ] ) > mytt ) {
     MSG( 0 , "Dt = +INF" << std::endl );
     return( HpINF );
     }
    }
   break;
  case quadratic:
   for( Index n = 0 ; n < CrrSGLen ; n++ )
    value += csol[ dict_stab[ n ] ] * csol[ dict_stab[ n ] ];
   value *= ( ( tt == t ? tt : ( t * t ) / tt ) / 2 );
   MSG( 0 , "Dt = " << value << std::endl );
   break;
  default:
   throw( NDOException( "OSIMPSolver::ReadDt: undecided stabilization" ) );
   }

 MSG( 0 , "Dt = " << value << std::endl );
 return( value );

 }  // end( OSIMPSolver::ReadDt )

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadSigma( cIndex wFi )
{
 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif
 const double *obj = osiSlvr->getObjCoefficients();

 // compute sigma  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 HpNum value = 0;
 if( wFi > NrFi ) {  // aggregate the linearization error and the rhs of the
                     // constraints of the difficult components

  for( Index i = 0 ; i < item_maxname ; i++ )
   if( ( dict_item[ i ] < InINF ) && WasInMP( i ) )
    value -= csol[ dict_item[ i ] ] * obj[ dict_item[ i ] ];

  // add the easy components cost part and the lower bounds contribution of
  // the difficult components  - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // note that is not necessary to divide any part of Sigma by rho

  for( Index i = 0 ; ++i <= NrFi ; )
   if( weasy[ i ] )
    for( Index j = 0; j < FIO->GetBNC( i ) ; j++ )
     value -= csol[ comp_col[ i ] + j ] * obj[ comp_col[ i ] + j ];
   else
    value -= csol[ comp_col[ i ] ] * obj[ comp_col[ i ] ];

  // If a global LB is present, we have LB * r in the objective function; but
  // we don't actually have r in the problem, rather rho = 1 - r. Therefore
  // we have LB * r = LB ( 1 - rho ) = - LB * rho + LB, and that's why in
  // the objective function we have - LB. But this also means that
  // *a fixed term LB has to be added to Sigma* if the LB is present

  if( LwrBnds[ 0 ] > -HpINF )
   value += obj[ comp_col[ 0 ] ] * ( 1 - csol[ comp_col[ 0 ] ] );

  // finally add the contribution of the slack variables - - - - - - - - - -
  // - - - - - -  - - - - - - -  - - - - - - - - - - - - - - - - - - - - - -

  if( stab != boxstep )
   for( Index i = 0; i < CrrSGLen ; i++ ) {
    int offset = 0;  // if the lower bound exists, there is a slack s_i^+
    if( Lower[ i ] > -LMINF ) {
     value -= Lower[ i ] * csol[ dict_slack[ i ] ];
     offset = 1;
     }

    if( Upper[ i ] < LMINF )
     value += Upper[ i ] * csol[ dict_slack[ i ] + offset ];
    }
  else
   for( Index i = 0; i < CrrSGLen ; i++ ) {
    // if -t < Lower[ i ], the stabilization constraint -t <= d[ i ] is
    // redundant, and the opposite for Upper[ i ]
    if( -t < Lower[ i ] )
     value -= Lower[ i ] * csol[ dict_slack[ i ] ];
    if( t > Upper[ i ] )
     value += Upper[ i ] * csol[ dict_slack[ i ] + 1 ];
    }

  MSG( 0 , "Sigma* + Sigma_L = " << value << std::endl );
  }
 else {
  if( weasy[ wFi ] )
   for( Index j = 0; j < FIO->GetBNC( wFi ) ; j++ )
    value -= csol[ comp_col[ wFi ] + j ] * obj[ comp_col[ wFi ] + j ];
  else {
   for( Index i = 0 ; i < item_maxname ; i++ )
    if( ( WComponent( i ) == wFi ) && WasInMP( i ) )
     value -= csol[ dict_item[ i ] ] * obj[ dict_item[ i ] ];
   value -= csol[ comp_col[ wFi ] ] * obj[ comp_col[ wFi ] ];
   }

  MSG( 0 , "Sigma*(" << wFi << ") = " << value << std::endl );
  }

 return( value );

 }  // end( OSIMPSolver::ReadSigma )

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadDStart( cHpNum tt )
{
 MSG( 0 , "OSIMPSolver::ReadDStart()" << std::endl );

 HpNum total = 0;
 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 switch( stab ) {
  case none: break;  // note: we assume the solution to be feasible
  case boxstep:
   for( Index i = 0 ; i < CrrSGLen ; i++ ) {
    if( -t >= Lower[ i ] )
     total += tt * csol[ dict_slack[ i ] ];

    if( t <= Upper[ i ] )
     total += tt * csol[ dict_slack[ i ] + 1 ];
    }
   break;
  case quadratic:
   for( Index i = 0 ; i < CrrSGLen ; i++ )
    total += csol[ dict_stab[ i ] ] * csol[ dict_stab[ i ] ] ;
   total *= ( tt / 2 );
   break;
  default:
   throw( NDOException( "OSIMPSolver::ReadDStart: undecided stabilization" )
	  );
  }

 MSG( 0 , "D* = " << total << std::endl );
 return( total );

 }  // end( OSIMPSolver::ReadDStart )

/*--------------------------------------------------------------------------*/

cLMRow OSIMPSolver::Readd( bool Fulld )
{
 MSG( 0 , "OSIMPSolver::Readd()\n" );

 if( Fulld ) {
  #if( PRESERVE_OSI_SOLS == 0 )
   const double *rsol = osiSlvr->getRowPrice();
  #endif

  return( rsol + comp_row[ NrFi ] );
  }
 else
  throw( NDOException( "OSIMPSolver::Readd: Fulld( false )" ) );

 }  // end( OSIMPSolver::Readd )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ReadZ( LMRow tz , cIndex_Set &I , Index &D , cIndex wFi )
{
 MSG( 0 , "OSIMPSolver::ReadZ() \n");

 // give it out in "dense" format  - - - - - - - - - - - - - - - - - - - - -

 D = CrrSGLen;
 I = 0;

 if( wFi == 0 ) { // return the gradient of 0-th component - - - - - - - - -
  VectAssign( tz , LMRow( 0 ) , CrrSGLen );
  for( Index i = RhoColSgPos ; i < RhoColBDm  ; i++ )
   tz[ RhoColBse[ i ] - comp_row[ NrFi ] ] = -RhoCol[ i ];
  return;
  }

 // make a convex combination of subgradient  - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 // use CoinPackedVector to simplify the work - - - - - - - - - - - - - - - -

 CoinPackedVector temp;
 if( wFi == 0 )    // return the gradient of 0-th component - - - - - - - - -
  temp = ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ 0 ] );
 else
  if( wFi > NrFi ) {  // take in account all the subgradients - - - - - - - -

   for( Index i = 0 ; i <= item_maxname ; i++ )
    if( ( IsSubG( i ) ) && WasInMP( i ) )
     temp = temp + csol[ dict_item[ i ] ] *
                  ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

   // add the subgradient of the easy components- - - - - - - - - - - - - - -

   for( Index i = 1 ; i <= NrFi ; i++ )
    if( weasy[ i ] )
     for( int j = comp_col[ i ] ; j < comp_col[ i ] + FIO->GetBNC( i ) ; j++ )
      temp = temp + csol[ j ] * ( osiSlvr->getMatrixByCol() )->getVector( j );

   if( wFi == InINF )   // add the linear part
    temp = temp + ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ 0 ] );
   }
  else { // disaggregated case  - - - - - - - - - - - - - - - - - - - - - - -
   if( weasy[ wFi ] )
    throw( NDOException( "OSIMPSolver::ReadZ: calling for easy component" ) );

   for( Index i = 0 ; i < item_maxname ; ++i )
    if( IsSubG( i ) && ( WComp( i ) == wFi ) && WasInMP( i ) )
     temp = temp + csol[ dict_item[ i ] ] *
      ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );
   }

 // assign tz the subgradient  - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 VectAssign( tz , LMNum( 0 ) , CrrSGLen );

 for( int j = 0 ; j < temp.getNumElements() ; j++ )
  if( temp.getIndices()[ j ] >= int( comp_row[ NrFi ] ) )
   tz[ temp.getIndices()[ j ] - comp_row[ NrFi ] ] = -temp.getElements()[ j ];

 }  // end( OSIMPSolver::ReadZ )

/*--------------------------------------------------------------------------*/

cHpRow OSIMPSolver::ReadMult( cIndex_Set &I , Index &D , cIndex wFi ,
			      const bool IncldCnst )
{
 MSG( 0 , "OSIMPSolver::ReadMult()\n" );

 if( wFi == 0 )
  throw( NDOException( "OSIMPSolver::ReadMult for 0-th component" ) );

 // Return the multipliers of the items of the difficult part. Don't
 // divide them by rho. In fact, the following holds:
 // \f$\sum_{k} \theta_k^i + \gamma^i + ( 1 - \rho ) = 1 \f$
 // where \f$ \gamma^i \f$ and \f$ ( 1 - \rho ) \f$ are, respectively,
 // the multiplier of the lower bound of the component wFi and
 // the global lower bound multiplier.

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 HpNum rho = csol[ comp_col[ 0 ] ];

 if( rho == 0 ) {  // in the case rho = 0, return a null vector- - - - - - -
  resizeI( 1 );
  tempI[ D = 0 ] = InINF;
  I = reinterpret_cast< cIndex_Set >( tempI );
  return( tempHP );
  }

 if( wFi <= NrFi )
  if( weasy[ wFi ] ) {  // easy component part- - - - - - - - - - - - - - -
   I = 0;
   D = FIO->GetBNC( wFi );
   resizeHP( D );
   VectAssign( tempHP , csol + comp_col[ wFi ] , D );
   }
  else {                // difficult component part - - - - - - - - - - - -
   resizeHP( MaxBSize );
   resizeI( MaxBSize + 1 );
   I = reinterpret_cast< cIndex_Set >( tempI );
   D = 0;
   if( IncldCnst ) {
    for( Index name = 0 ; name < item_maxname ; name++ )
     if( ( WComponent( name ) == wFi ) && WasInMP( name ) )
      if( ( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 )
       tempI[ D++ ] = name;
    }
   else
    for( Index name = 0 ; name < item_maxname ; name++ )
     if( IsSubG( name ) && ( WComp( name ) == wFi ) && WasInMP( name ) )
      if( ( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 )
       tempI[ D++ ] = name;

   tempI[ D ] = InINF;
   }
 else {  // wFi == InINF- - - - - - - - - - - - - - - - - - - - - - - - - -
  resizeHP( MaxBSize );
  resizeI( MaxBSize + 1 );
  I = reinterpret_cast< cIndex_Set >( tempI );
  D = 0;
  if( IncldCnst ) {
   for( Index name = 0 ; name < item_maxname ; name++ )
    if( ( dict_item[ name ] < InINF ) && WasInMP( name ) )
     if( ( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 )
      tempI[ D++ ] = name;
   }
  else
   for( Index name = 0 ; name < item_maxname ; name++ )
    if( IsSubG( name ) && WasInMP( name ) )
     if( ( tempHP[ D ] = csol[ dict_item[ name ] ] ) != 0 )
      tempI[ D++ ] = name;

  tempI[ D ] = InINF;
  }

 return( tempHP );

 }  // end( OSIMPSolver::ReadMult )

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadLBMult( cIndex wFi )
{
 MSG( 0 , "OSIMPSolver::ReadLBMult()\n" );

 if( wFi == 0 )
  throw( NDOException( "OSIMPSolver::ReadLBMult for 0-th component" ) );

 // return the multipliers of the lower bounds - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *csol = osiSlvr->getColSolution();
 #endif

 HpNum value;
 if( wFi <= NrFi ) {  // return \f$ \gamma^i \f$ - - - - - - - - - - - - - - -
  if( weasy[ wFi ] )  // - - - - - - - - - - - - - - - - - - - - - - - - - - -
   throw( NDOException( "OSIMPSolver::ReadLBMult for an easy component" ) );

  value = csol[ comp_col[ wFi ] ];
  }
 else  // return  \f$r = ( 1 - \rho ) \f$  - -  - - - - - - - - - - - - - - -
  value = 1 - csol[ comp_col[ 0 ] ];

 return( value );
 }

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadGid( cIndex Nm )
{
 MSG( 0 , "OSIMPSolver::ReadGid() called\n" );

 HpNum value = 0;
 if( Nm >= MaxBSize ) {  // 0-th component: return delta * b
  #if( PRESERVE_OSI_SOLS == 0 )
   const double *rsol = osiSlvr->getRowPrice();
  #endif
  for( Index i = RhoColSgPos ; i < RhoColBDm ; i++ )
   value -= rsol[ RhoColBse[ i ] ] * RhoCol[ i ];
  MSG( 0 , "b * delta = " << value << std::endl );
  }
 else {
  if( dict_item[ Nm ] == InINF )
   throw( NDOException( "OSIMPSolver::ReadGid: unused item name" ) );

  // note: the scalar product is available for all items, even those that
  //       were not present in the last master problem (i.e., such that
  //       wcomp[ Nm ] & LLBIndex == true) because it is computed in
  //       Check[SubG/Const]() and saved in SetItem().

  value = GiPerd[ Nm ];

  MSG( 0 , "Item[" << Nm << "] * d = " << value << std::endl );
  }

 return( value );

 }  // end( OSIMPSolver::ReadGid )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::MakeLambda1( cHpRow Lmbd , HpRow Lmbd1 , cHpNum Tau )
{
 MSG( 0, "OSIMPSolver::MakeLambda1()\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *primal = rsol + comp_row[ NrFi ];

 if( useactiveset )   // Lambda and Lambda1 are in sparse format- - - - - -
  for( Index i = 0 ; i < Asetdim ; i++ ) {
   HpNum step = ( Tau / t ) * primal[ Aset[ i ] ];
   if( step > Upper[ i ] ) step = Upper[ Aset[ i ] ];
   if( step < Lower[ i ] ) step = Lower[ Aset[ i ] ];
   Lmbd1[ i ] = Lmbd[ i ] + step;
   }
 else   // Lambda and Lambda1 are in dense format- - - - - - - - - - - - - -
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   HpNum step = ( Tau / t ) * primal[ i ];
   if( step > Upper[ i ] ) step = Upper[ i ];
   if( step < Lower[ i ] ) step = Lower[ i ];
   Lmbd1[ i ] = Lmbd[ i ] + step;
   }

 }  // end( OSIMPSolver::MakeLambda1 )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SensitAnals( HpNum &lp , HpNum &cp )
{
 throw( NDOException( "OSIMPSolver::SensitAnals not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE PROBLEM ----------------*/
/*--------------------------------------------------------------------------*/

Index OSIMPSolver::BSize( cIndex wFi )
{
 if( ( wFi > 0 ) && ( wFi <= NrFi ) )
  return( NSubG[ wFi ] + NConst[ wFi ] );
 else
  return( NSubG[ 0 ] + NConst[ 0 ] );

 }  // end( OSIMPSolver::BSize )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::BCSize( cIndex wFi )
{
 if( ( wFi > 0 ) && ( wFi <= NrFi ) )
  return( NConst[ wFi ] );
 else
  return( NConst[ 0 ] );

 }  // end( OSIMPSolver::BCSize )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::MaxName( cIndex wFi )
{
 Index i = item_maxname;
 if( wFi != InINF )
  for( ; i > 0 ; i-- )
   if( WComponent( i - 1 ) == wFi )
    break;

 return( i );

 }  // end( OSIMPSolver::MaxName )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::WComponent( cIndex i )
{
 return( wcomp[ i ] < InINF ? WComp( i ) : InINF );

 }  // end( OSIMPSolver::WComponent )

/*--------------------------------------------------------------------------*/

bool OSIMPSolver::IsSubG( cIndex i )
{
 return( wcomp[ i ] < InINF ? bool( wcomp[ i ] & LBIndex ) : false );

 }  // end( OSIMPSolver::IsSubG )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::NumNNVars( void )
{
 return( NNVars );

 }  // end( OSIMPSolver::NumNNVars )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::NumBxdVars( void )
{
 return( BxdVars );

 }  // end( OSIMPSolver::NumBxdVars )

/*--------------------------------------------------------------------------*/

bool OSIMPSolver::IsNN( cIndex i )
{
 return( Lower[ i ] > -LMINF );

 }  // end( OSIMPSolver::IsNN )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::CheckIdentical( const bool Chk )
{
 checkID = Chk;

 }  // end( OSIMPSolver::CheckIdentical )

/*--------------------------------------------------------------------------*/

cHpRow OSIMPSolver::ReadLinErr( void )
{
 resizeHP( MaxBSize );
 VectAssign( tempHP , HpNum( 0 ) , MaxBSize );

 const double *linerr = osiSlvr->getObjCoefficients();
 for( Index i = 0 ; i < MaxBSize ; i++ )
  if( dict_item[ i ] < InINF )
   tempHP[ i ] = - linerr[ dict_item[ i ] ];  // change sign!

 return( tempHP );

 }  // end( OSIMPSolver::ReadLinErr )

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::ReadLowerBound( cIndex wFi )
{
 if( ( wFi <= NrFi ) && weasy[ wFi ] )
  throw( NDOException( "OSIMPSolver::ReadLowerBound: wFi easy component" ) );

 return( LwrBnds[ wFi > NrFi ? 0 : wFi ] );
 }

/*--------------------------------------------------------------------------*/

HpNum OSIMPSolver::EpsilonD( void )
{
 return( FsbEps );

 }  // end( OSIMPSolver::EpsilonD )

/*--------------------------------------------------------------------------*/

bool OSIMPSolver::FiBLambdaIsExact( cIndex wFi )
{
 return( weasy[ wFi ] );

 }  // end( OSIMPSolver::FiBLambdaIsExact )

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR ADDING / REMOVING / CHANGING DATA ---------------*/
/*--------------------------------------------------------------------------*/

SgRow OSIMPSolver::GetItem( cIndex wFi )
{
 NewItemFi = wFi;  // mark the name of the component
 return( NewItem );

 }  // end( OSIMPSolver::GetItem )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetItemBse( cIndex_Set SGBse , cIndex SGBDm )
{
 if( ( NewItemBse = SGBse ) )
  NewItemBDm = SGBDm;    // sparse format
 else
  NewItemBDm = CrrSGLen; // dense format

 }  // end( OSIMPSolver::SetItemBse )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::CheckSubG( cHpNum DFi , cHpNum Tau , HpNum &Ai ,
			      HpNum &ScPri )
{
 MSG( 0 , "OSIMPSolver::CheckSubG() called\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *primal = rsol + comp_row[ NrFi ];

 // check if the subgradient is identical to any other in the bundle - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index IsIde = CheckBCopy();

 // compute the scalar product Gi^{top} delta- - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Tau == 0 )     // subgradient in Lamda, ScPri not needed
  ScPri = 0;
 else
  if( NewItemBse )  // sparse format
   ScPri = ScalarProduct( NewItem , primal , NewItemBse );
  else              // dense format
   ScPri = ScalarProduct( NewItem , primal , CrrSGLen );

 // compute the linearization error of the new subgradient - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // This is done only if Tau > 0, as Tau <= 0 signals that Ai is already the
 // linearization error in the current point

 if( Tau > 0 )
  Ai = Ai - ( DFi - ( Tau / t ) * ScPri );

 NewItemprice = Ai;
 NewItemScPri = ScPri;
 NewItemisSG = true;   // this item is a subgradient ...

 return( IsIde );

 } // end( OSIMPSolver::CheckSubG )

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::CheckCnst( HpNum &Ai , HpNum &ScPri , cHpRow CrrPnt )
{
 MSG( 0 , "OSIMPSolver::CheckCnst() called\n" );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *primal = rsol + comp_row[ NrFi ];

 // check if the constraint is identical to any other in the bundle- - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index IsIde = CheckBCopy();

 // compute ScPri and right hand side of the constraint- - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( NewItemBse ) {  // sparse format- - - - - - - - - - - - - - - - - - -
  ScPri = ScalarProduct( NewItem , primal , NewItemBse );

  // compute the right hand side of the constraint
  if( Aset )         // newitem is sparse & is using the active set
   Ai -= ScalarProduct( NewItem , NewItemBse , CrrPnt , Aset );
  else               // newitem is sparse & is not using the active set
   Ai -= ScalarProduct( NewItem , CrrPnt , NewItemBse );
  }
 else {              // dense format - - - - - - - - - - - - - - - - - - -
  ScPri = ScalarProduct( NewItem , primal , CrrSGLen );

  // compute the right hand side of the constraint
  if( Aset )         // newitem is dense & activeset
   Ai -= ScalarProduct( CrrPnt , NewItem , Aset );
  else               // newitem is dense & not activeset
   Ai -= ScalarProduct( NewItem , CrrPnt , CrrSGLen );
  }

 NewItemprice = Ai;
 NewItemScPri = ScPri;
 NewItemisSG = false;  // this item is a constraint ...

 return( IsIde );

 }  // end( OSIMPSolver::CheckCnst )

/*--------------------------------------------------------------------------*/

bool OSIMPSolver::ChangesMPSol( void )
{
 MSG( 0 , "OSIMPSolver::ChangesMPSol() called\n" );

 if(! NewItemFi )
  throw( NDOException(
           "OSIMPSolver::ChangesMPSol: a calling to 0-th is not allowed" ) );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 bool ChgMP = true;

 if( NewItemisSG ) {
  /* check whether or not the new subgradient provides a cut in Master
     Problem, that is verify v >= G1 * dir - Alfa1 is violated */

  if( NewItemFi <= NrFi ) {
   if( weasy[ NewItemFi ] )
    throw( NDOException(
     "OSIMPSolver::ChangesMPSol: check for an easy comp. is not allowed" )
     );
   }
  else
   throw( NDOException( "OSIMPSolver::ChangesMPSol: "
		  "check for an aggregated subgradient is not allowed" ) );

  if( rsol[ comp_row[ NewItemFi - 1 ] ] + NewItemprice - NewItemScPri >=
      - FsbEps * std::max( ABS( NewItemprice ) , double(1) ) )
   ChgMP = false;

  }
 else // check whether or not the new constraint provides a cut
      // in Master Problem, that is verify Alfa1 >= G1 * dir is violated
  if( NewItemprice - NewItemScPri >=
      - FsbEps * std::max( ABS( NewItemprice ) , double( 1 ) ) )
   ChgMP = false;

 return( ChgMP );

 } // end( OSIMPSolver::ChangesMPSol )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetItem( cIndex Nm )
{
 // note: the columns of the coefficient matrix (for the relevant part)
 //       contain the *opposite* of the subgradient

 MSG( 0 , "OSIMPSolver::SetItem() \n");

 if( ( ! NewItemFi ) && ( Nm == InINF ) ) {  // the 0th component
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // we have to update the rho column, in particular the part corresponding
  // to the subgradient (while leaving the rest untouched)
  RhoColBDm = RhoColSgPos;

  if( NewItemBse )  // sparse format - - - - - - - - - - - - - - - - - - - -
   for( Index i = 0 ; i < NewItemBDm ; i++ ) {
    if( NewItem[ i ] ) {
     RhoColBse[ RhoColBDm ] = comp_row[ NrFi ] + NewItemBse[ i ];
     RhoCol[ RhoColBDm++ ] = - NewItem[ i ];
     }
    }
  else 	            // dense format  - - - - - - - - - - - - - - - - - - - -
   for( Index i = 0 ; i < CrrSGLen ; i++ )
    if( NewItem[ i ] ) {
     RhoColBse[ RhoColBDm ] = comp_row[ NrFi ] + i;
     RhoCol[ RhoColBDm++ ] = - NewItem[ i ];
     }

  UpdateRhoCol();  // remove the previous version and add the new one
  }
 else {  // it's an actual item- - - - - - - - - - - - - - - - - - - - - - -

  resizeI( NewItemBDm + 1 );   // note: NewItemBD == CrrSGLen in the dense
  resizeHP( NewItemBDm + 1 );  // case

  // copy the item in the temporary structure  - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index i = 0;
  if( NewItemBse )
   for( ; i < NewItemBDm ; i++ ) {
    tempI[ i ] = NewItemBse[ i ] + comp_row[ NrFi ];
    tempHP[ i ] = - NewItem[ i ];
    }
  else
   for( ; i < CrrSGLen ; i++ ) {
    tempI[ i ] = i + comp_row[ NrFi ];
    tempHP[ i ] = - NewItem[ i ];
    }

  if( NewItemisSG ) {  // the new item is a subgradient

   // if the new subgradient is the first one for its component, the
   // convexity constraint for that component has to be activated, unless an
   // individual lower bound is set, in which case \gamma_i must be left free
   // if \gamma_i is fixed to 0, also set its cost to 0; if it is not fixed
   // its cost, which depends on whether or not a global lower bound has also
   // been set, does not change
   if( ( ! NSubG[ NewItemFi ] ) && ( ( LwrBnds[ NewItemFi ] == -HpINF ) ) ) {
    osiSlvr->setColUpper( comp_col[ NewItemFi ] , 0 );
    osiSlvr->setObjCoeff( comp_col[ NewItemFi ] , 0 );
    }

   // add the multiplier to the simplex set   
   tempI[ i ] = comp_row[ NewItemFi - 1 ];
   tempHP[ i++ ] = 1;
   NSubG[ NewItemFi ]++;
   (*NSubG)++;

   // set name, plus it's a "new" item, plus is a subgradient
   wcomp[ Nm ] = NewItemFi | IMask;
   }
  else {               // the new item is a constraint
   NConst[ NewItemFi ]++;
   (*NConst)++;

   wcomp[ Nm ] = NewItemFi | LLBIndex;  // set name, plus it's a "new" item
   }

  // also record the scalar product for future use (ChangeCurrPoint) - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  GiPerd[ Nm ] = NewItemScPri;

  // put the item in the dictionary  - - - - - - - - - - - - - - - - - - - -

  dict_item[ Nm ] = osiSlvr->getNumCols();

  // add the column to the master problem- - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->addCol( int( i ) , tempI , tempHP , 0 , osiSlvr->getInfinity() ,
		   - NewItemprice );

  #if( PRESERVE_OSI_SOLS )
   // ensure that csol[] and rcost[] are of the right size - - - - - - - - -
   if( dict_item[ Nm ] >= csols ) {
    // allocate new vectors a bit larger to avoid doing this too often
    auto ncsol = new double[ dict_item[ Nm ] + 20 ];
    auto nrcst = new double[ dict_item[ Nm ] + 20 ];

    // copy over previous solution and initialize new elements somehow
    VectAssign( ncsol , csol , csols );
    VectAssign( nrcst , rcst , csols );
    ncsol[ dict_item[ Nm ] ] = nrcst[ dict_item[ Nm ] ] = 0;

    // delete previous versions
    delete[] csol;
    delete[] rcst;

    // update size
    csols = dict_item[ Nm ] + 20;

    // record new pointers
    csol = ncsol;
    rcst = nrcst;
    }
  #endif

  MSG( 0 , "Item " << Nm << " of the component " << NewItemFi
	<< " has the entry " << dict_item[ Nm ] << std::endl);

  if( item_maxname < Nm + 1 )
   item_maxname = Nm + 1;
  }
 }  // end( OSIMPSolver::SetItem )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SubstItem( cIndex Nm )
{
 MSG( 0 , "OSIMPSolver::SubstItem()\n" );

 //if( - osiSlvr->getObjCoefficients()[ dict_item[ Nm ] ] > NewItemprice )
 osiSlvr->setObjCoeff( dict_item[ Nm ] , - NewItemprice );

 } // end( OSIMPSolver::SubstItem )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::RmvItem( cIndex i )
{
 const int index = dict_item[ i ];  // mark the column deleted - - - - - - -
 osiSlvr->deleteCols( 1 , &index );

 #if( PRESERVE_OSI_SOLS )
  // having deleted the column 'index' from the osiSlvr, all the stored
  // solution information corresponding to columns (csol[] and rcst[]) must
  // be shifted left by one from position index onwards

  ShiftVect( csol + index , csols - index - 1 );
  ShiftVect( rcst + index , csols - index - 1 );
 #endif

 Index wFi = WComp( i );

 if( IsSubG( i ) ) {
  NSubG[ wFi ]--;
  (*NSubG)--;
  }
 else {
  NConst[ wFi ]--;
  (*NConst)--;
  }

 // update the item's vocabulary - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 dict_item[ i ] = wcomp[ i ] = InINF;

 while( item_maxname && ( dict_item[ item_maxname - 1 ] == InINF ) )
  --item_maxname;

 // if the deleted subgradient is the last one for its component, the
 // convexity constraint for that component has to be deactivated; however,
 // if an individual bound is set, this is not necessary since \gamma_i is
 // free already (and it would be wrong to set its cost to 0!)
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ! NSubG[ wFi ] ) && ( LwrBnds[ wFi ] == -HpINF ) ) {
  osiSlvr->setColUpper( comp_col[ wFi ] , osiSlvr->getInfinity() );

  // the cost of \gamma_{wFi} is 0 if no global lower bound is set, otherwise
  // is set < than (-) that of \rho
  osiSlvr->setObjCoeff( comp_col[ wFi ] ,
			( LwrBnds[ 0 ] == -HpINF ) ? 0 : LwrBnds[ 0 ] - 1 );
  }

 // update the dictionaries- - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index Nm = 0 ; Nm < item_maxname ; Nm++ )
  if( ( dict_item[ Nm ] != InINF ) &&
      ( dict_item[ Nm ] > Index( index ) ) )
   dict_item[ Nm ]--;

 for( Index Nm = 0 ; Nm < CrrSGLen ; Nm++ ) {
  if( ( dict_slack[ Nm ] != InINF ) &&
      ( dict_slack[ Nm ] > Index( index ) ) )
   dict_slack[ Nm ]--;
  if( stab == quadratic )
   if( ( dict_stab[ Nm ] != InINF ) &&
       ( dict_stab[ Nm ] > Index( index ) ) )
    dict_stab[ Nm ]--;
  }

 if( comp_col[ 0 ] > Index( index ) )
  comp_col[ 0 ]--;

 MSG( 0 , "removed item " << i << " of component " << wFi << std::endl );

 }  // end( OSIMPSolver::RmvItem )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::RmvItems( void )
{
 MSG( 0 , "OSIMPSolver::RmvItems()\n" );

 if( ! item_maxname )  // there is nothing to do ...
  return;

 VectAssign( NSubG , Index( 0  ) , NrFi + 1 );
 VectAssign( NConst , Index( 0 ) , NrFi + 1 );

 int count = 0;
 resizeI( item_maxname );

 for( Index Nm = 0 ; Nm < item_maxname ; Nm++ )
  if( dict_item[ Nm ] != InINF )
   tempI[ count++ ] = dict_item[ Nm ];

 // delete the columns of all the items- - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 osiSlvr->deleteCols( count , tempI );

 // fill the vectors for the shifting of the vocabularies  - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index elem;
 for( Index Nm = 0 ; Nm < item_maxname ; Nm++ )
  if( dict_item[ Nm ] != InINF ) {
   elem = dict_item[ Nm ];

   // update the left part of tempI  - - - - - - - - - - - - - - - - - - - -

   for( Index i = Nm + 1; i < item_maxname ; i++ )
    if( ( dict_item[ i ] > elem ) && ( dict_item[ i ] != InINF ) )
     dict_item[ i ]--;

   // update the position of rho column  - - - - - - - - - - - - - - - - - -

   if( comp_col[ 0 ] > elem )
    comp_col[ 0 ]--;

   // update the slack and stabilization dictionaries  - - - - - - - - - - -

   for( Index i = 0 ; i < CrrSGLen ; i++ ) {
    if( ( dict_slack[ i ] != InINF ) && ( dict_slack[ i ] > elem ) )
     dict_slack[ i ]--;
    if( stab == quadratic )
     if( ( dict_stab[ i ] != InINF ) && ( dict_stab[ i ] > elem ) )
      dict_stab[ i ]--;
    }

   } // end check  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // for each component that does *not* have an individual lower bound set,
 // deactivate the convexity constraint; if the individual lower bound is set
 // this is not necessary since \gamma_i is free already (and it would be
 // wrong to set its cost to 0!)
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i++ < NrFi ; )
  if( LwrBnds[ i ] == -HpINF ) {
   osiSlvr->setColUpper( comp_col[ i ] , osiSlvr->getInfinity() );

   // the cost of \gamma_i is 0 if no global lower bound is set, otherwise
   // is set < than (-) that of \rho
   osiSlvr->setObjCoeff( comp_col[ i ] ,
			 ( LwrBnds[ 0 ] == -HpINF ) ? 0 : LwrBnds[ 0 ] - 1 );
   }

 VectAssign( dict_item , InINF , item_maxname );
 VectAssign( wcomp , InINF , item_maxname );
 item_maxname = 0;

 }  // end( OSIMPSolver::RmvItems )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::SetActvSt( cIndex_Set AVrs , cIndex AVDm )
{
 MSG( 0 , "OSIMPSolver::SetActvSt()\n" );

 if( ! useactiveset )
  throw( NDOException( "OSIMPSolver::SetActvSt(): Active set not declared"
		       ) );

 if( ! AVrs ) {	 //  all the variables will be not active- - - - - - - - - -
  while( Asetdim > 0 )
   deactivate( Aset[ Asetdim-- ] );
  }
 else
  if( ! Aset ) {  // if Aset is empty, set all the variables in AVrs as active
    while ( Asetdim < AVDm )
     activate( AVrs[ Asetdim++ ] );
   }
  else
   for( cIndex *newA = AVrs ; ( *newA < InINF ) ||
	                      ( *Aset < InINF ) ; ) {
    if( *newA < *Aset )	{
     activate( *newA );
     newA++;
     }
    else
     if( *newA > *Aset ) {
      deactivate( *Aset );
      Aset++;
      }
     else {
      newA++;
      Aset++;
      }
    }  // end( for )

 Aset = AVrs;
 Asetdim = AVDm;

 }  // end( OSIMPSolver::SetActvSt )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::AddActvSt( cIndex_Set Addd , cIndex AdDm , cIndex_Set AVrs )
{
 MSG( 0 , "OSIMPSolver::AddActvSt() called\n" );

 if( ! useactiveset )
  throw( NDOException( "OSIMPSolver::AddActvSt(): Active set not declared"
		       ) );
 Aset = AVrs;
 Asetdim += AdDm;

 for( ; *Addd < InINF ; Addd++ )
  activate( *Addd );

 }  // end( OSIMPSolver::AddActvSt )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::RmvActvSt( cIndex_Set Rmvd , cIndex RmDm , cIndex_Set AVrs )
{
 MSG( 0 , "OSIMPSolver::RmvActvSt() called\n" );

 Aset = AVrs;
 Asetdim = Asetdim - RmDm;

 for( ; *Rmvd < InINF ; Rmvd++ )
  deactivate( *Rmvd );

 } // end( OSIMPSolver::RmvActvSt )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::AddVars( cIndex NNwVrs )
{
 MSG( 0 , "OSIMPSolver::AddVars() called\n" );

 if( ! NNwVrs )  // adding 0 new variables
  return;        // all is done

 cIndex NewSGLen = CrrSGLen + NNwVrs;
 if( NewSGLen > MaxSGLen )
  throw( NDOException( "OSIMPSolver::AddVars: too many variables" ) );

 // allocate the memory for new rows - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index size = 0;
 int * rowStarts = new int[ NNwVrs + 1 ];      // start point of the rows
 HpRow rowlb = new HpNum[ NNwVrs ];            // lower bound of constraints
 HpRow rowub = new HpNum[ NNwVrs ];            // upper bound of constraints

 // recovery the description of  matrix A[ wFi ] of the "easy" - - - - - - -
 // components - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SIndex_Mat Abeg = new SIndex_Set[ NrFi ];
 SIndex_Mat Aind = new SIndex_Set[ NrFi ];
 Mat Aval = new Row[ NrFi ];

 Index_Set BNC = new Index[ NrFi ];
 Index_Set ANZ = new Index[ NrFi ];

 for( Index j = 1 ; j <= NrFi ; j++)
  if( ( BNC[ j - 1 ] = FIO->GetBNC( j ) ) ) {
   ANZ[ j - 1 ] = FIO->GetANZ( j , CrrSGLen , NewSGLen );
   if( ANZ[ j - 1 ] ) {
	size += ANZ[ j - 1 ];
    Abeg[ j - 1 ] = new SIndex[ BNC[ j - 1 ] + 1 ];
    Aind[ j - 1 ] = new SIndex[ ANZ[ j - 1 ] ];
    Aval[ j - 1 ] = new Number[ ANZ[ j - 1 ] ];
    FIO->GetADesc( j , Abeg[ j - 1 ] , Aind[ j - 1 ] , Aval[ j - 1 ] ,
		   CrrSGLen , NewSGLen );
    }
   }
  else {
   Abeg[ j - 1 ] = Aind[ j - 1 ] = 0;
   Aval[ j - 1 ] = 0;
   ANZ[ j - 1 ] = 0;
   }

 // recovery all part of sub-gradients corresponding to variables in Lambda
 // from index CrrSGLen to ( NewSGLen - 1 ) - - - - - - - - - - - - - - - - -

 cIndex_Set SGBse;
 SgMat tempG = new SgRow[ item_maxname ];

 for( Index j = 0 ; j < item_maxname ; j++ )
  if( dict_item[ j ] < InINF ) {  // ask the components [ CrrSGLen ,
   tempG[ j ] = new SgNum[ NNwVrs ];     // NewSGLen ) for item j
   cIndex SGBDim = FIO->GetGi( tempG[ j ] , SGBse , j , CrrSGLen , NewSGLen );
   if( SGBse )
    Densify( tempG[ j ] - CrrSGLen , SGBse , SGBDim , NewSGLen , CrrSGLen );

   size += NNwVrs;
   }
  else
   tempG[ j ] = 0;

 // recovery the subgradient of the linear 0-th component of Fi- - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 SgRow tempG0 = new SgNum[ NNwVrs ];
 cIndex SGBDim = FIO->GetGi( tempG0 , SGBse , FIO->GetMaxName() , CrrSGLen ,
			     NewSGLen );
 if( SGBse )
  Densify( tempG0 - CrrSGLen , SGBse , SGBDim , NewSGLen , CrrSGLen );

 size += NNwVrs;

 // construct the rows to the problem, one after another - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 resizeI( size );
 resizeHP( size );

 Index count = 0;
 for( Index i = CrrSGLen ; i < NewSGLen ; i++ ) {

  rowStarts[ i - CrrSGLen ] = count;

  // write the matrices A[ wFi ] in the problem, one for each easy component

  for( Index j = 1 ; j <= NrFi ; j++)
   if( ANZ[ j - 1 ] ) {
    for( Index k = 0 ; k < BNC[ j - 1 ] ; k++ )
     for( SIndex l =  Abeg[ j - 1 ][ k ]; l < Abeg[ j - 1 ][ k + 1 ] ; l++ )
      if( Index( Aind[ j - 1 ][ l ] ) == i ) {
       tempI[ count ] = k + comp_col[ j ];
       tempHP[ count++ ] = Aval[ j - 1 ][ l ];
       }
    }

  // write in both tempI and tempHP the subgradient part - - - - - - - - - -

  for( Index j = 0 ; j < item_maxname ; j++ )
   if( dict_item[ j ] < InINF ) {
    tempI[ count ] = dict_item[ j ];
    tempHP[ count++ ] = - tempG[ j ][ i - CrrSGLen ];
    }

  // write in both tempI and tempHP the subgradient of the linear 0-th
  // component

  RhoColBse[ RhoColBDm ] = comp_row[ NrFi ] + i;     // update the column
  RhoCol[ RhoColBDm++ ] = - tempG0[ i - CrrSGLen ];  // of rho

  tempI[ count ] = comp_col[ 0 ];
  tempHP[ count++ ] = -tempG0[ i - CrrSGLen ];
  }

 rowStarts[ NNwVrs ] = count;

 // add the rows to problem - - - - - - - - - - - - - - - - - - - - - - - -

 if( useactiveset ) {
  MSG( 0 , "Add constraints for the primal variables as inactive ones\n" );
  for( Index j = 0 ; j < NNwVrs ;  j++ ) {
   rowlb[ j ] = - osiSlvr->getInfinity();
   rowub[ j ] = osiSlvr->getInfinity();
   }
  osiSlvr->addRows( NNwVrs , rowStarts, tempI , tempHP , rowlb , rowub );
  }
 else {
  MSG( 0 , "Add constraints for the primal variables \n" );
  for( Index j = 0 ; j < NNwVrs ; ++j )
   rowlb[ j ] = rowub[ j ] = 0;

  osiSlvr->addRows( NNwVrs , rowStarts, tempI , tempHP , rowlb , rowub );
  }

 // deallocate the memory- - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] tempG0;
 for( Index j = 0 ; j < item_maxname ; j++ )
   delete[] tempG[ j ];
 delete[] tempG;

 delete[] BNC;
 delete[] ANZ;

 for( Index j = 1 ; j <= NrFi ; j++ ) {
  delete[] Abeg[ j - 1 ];
  delete[] Aind[ j - 1 ];
  delete[] Aval[ j - 1 ];
  }

 delete[] Abeg;
 delete[] Aind;
 delete[] Aval;

 delete[] rowStarts;
 delete[] rowlb;
 delete[] rowub;

 // add the slacks to the problem- - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index nc = osiSlvr->getNumCols();  // number of columns - 1
 HpRow obj;                         // coefficients of the objective
                                    // function relative to slacks
 int *columnStarts;
 if( dict_stab ) {
  columnStarts = new int[ ( 3 * NNwVrs ) + 1 ];
  resizeI( 3 * NNwVrs );
  resizeHP( 3 * NNwVrs );
  rowlb = new HpNum[ 3 * NNwVrs ];
  }
 else {
  columnStarts = new int[ ( 2 * NNwVrs ) + 1 ];
  resizeI( 2 * NNwVrs );
  resizeHP( 2 * NNwVrs );
  rowlb = new HpNum[ 2 * NNwVrs ];
  }

 // redefine the size both of tempHp and tempI - - - - - - - - - - - - - - -

 bool slack_p;
 bool slack_m;

 count = 0;
 for( Index i = CrrSGLen ; i < NewSGLen ; i++ ) {
  slack_p = false;
  slack_m = false;

  if( ( Upper[ i ] = FIO->GetUB( i ) ) != LMINF )
   slack_m = true;

  if( FIO->GetUC( i ) )
   Lower[ i ] = -LMINF;
  else {
   Lower[ i ] = 0;
   slack_p = true;
   }

  if( slack_p )
   NNVars++;

  if( slack_p || slack_m )
   BxdVars++;

  // on the basis of the adopted stabilization type we have
  // either one or two slacks  - - - - - - - - - - - - - - - - - - - - - - -

  switch ( stab ) {
   case none:    break;
   case boxstep: slack_p = slack_m = true;
                 break;
   case quadratic:
	columnStarts[ count ] = count;
	dict_stab[ i ] = count + nc;
	rowlb[ count ] = - osiSlvr->getInfinity();
	tempI[ count ] = comp_row[ NrFi ] + i;
	tempHP[ count++ ] = 1;
	MSG( 0 , "Add z_"<< i << std::endl );
	break;
   default:
    throw( NDOException( "OSIMPSolver::AddVars: undecided stabilization" ) );
   }

  if( slack_p ) {
   columnStarts[ count ] = count;
   dict_slack[ i ] = count + nc;
   rowlb[ count ] = 0;
   tempI[ count ] = comp_row[ NrFi ] + i;
   tempHP[ count++ ] = 1;
   MSG( 0 , "Add the slack s_"<< i <<"+\n" );
   }

  if( slack_m )	{
   columnStarts[ count ] = count;

   // if slack_p is present, slack_m is not indicated
   if( ! slack_p )
    dict_slack[ i ] = count + nc;

   rowlb[ count ] = 0;
   tempI[ count ] = comp_row[ NrFi ] + i;
   tempHP[ count++ ] = -1;
   MSG( 0 , "Add the slack s_" << i << "-\n" );
   }
  }

 columnStarts[ count ] = count;  // set the end both of tempHp and tempI
 rowub = new HpNum[ count ];
 obj = new HpNum[ count ];
 for( Index j = 0 ; j < count ; j++ ) {
  rowub[ j ] = osiSlvr->getInfinity();
  obj[ j ] = 0;
  }

 osiSlvr->addCols( count , columnStarts , tempI , tempHP , rowlb , rowub ,
		   obj );

 // deallocate the memory - - - - - - - - - - - - - - - - - - - - - - - - - -

 delete[] columnStarts;
 delete[] rowlb;
 delete[] rowub;
 delete[] obj;

 // update the current items length - - - - - - - - - - - - - - - - - - - - -

 CrrSGLen = NewSGLen;

 // change the price of the slack variables - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // change the prices for the stabilization
 tUpdatePrices( CrrSGLen - NNwVrs , CrrSGLen ); 
 if( stab != boxstep )  // change the slack prices, unless they have just
  ptUpdatePrices( CrrSGLen - NNwVrs , CrrSGLen );  // been changed 
 
 }  // end( OSIMPSolver::AddVars )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::RmvVars( cIndex_Set whch , Index hwmny )
{
 MSG( 0 , "OSIMPSolver::RmvVars() called\n" );

 int count;
 if( whch ) {  // removing a specific subset of variables - - - - - - - - - -
               // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Index newCrrSGLen = CrrSGLen - hwmny;

  if( dict_stab )
   resizeI( 3 * hwmny );
  else
   resizeI( 2 * hwmny );

  // remove the rows in the dynamic part - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < hwmny ; i++ )
   tempI[ i ] = whch[ i ] + comp_row[ NrFi ];

  osiSlvr->deleteRows( hwmny , tempI );

  // update rho column- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // note that the initial "static" part does not change, and therefore in
  // particular RhoColSgPos remains the same
  CoinPackedVector rhoCpy = 
                   ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ 0 ] );
  for( Index j = rhoCpy.getNumElements() ; j-- > 0 ; ) {
   RhoCol[ j ] = rhoCpy.getElements()[ j ];
   RhoColBse[ j ] = rhoCpy.getIndices()[ j ];
   }
  RhoColBDm = rhoCpy.getNumElements();

  // mark the columns to delete   - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  count = 0;
  for( Index i = 0 ; i < hwmny ; i++ ) {
   cIndex wi = whch[ i ];
   bool islwr =  ( Lower[ wi ] > -LMINF );
   bool isupr = ( Upper[ wi ] < LMINF );

   if( islwr )
    NNVars--;

   if( islwr || isupr )
    BxdVars--;

   if( dict_slack[ wi ] != InINF ) {
    tempI[ count++ ] = dict_slack[ wi ];
    if( ( islwr && isupr ) || ( stab == boxstep ) )
     tempI[ count++ ] = dict_slack[ wi ] + 1;
    }
   if( dict_stab )
    tempI[ count++ ] = dict_stab[ wi ];
   }

  // delete the columns of all the items- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->deleteCols( count , tempI );

  // shift the vocabularies - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Index elem;
  for( int Nm = 0 ; Nm < count ; Nm++ ) {
   elem = tempI[ Nm ];

   // update the left part of tempI  - - - - - - - - - - - - - - - - - - - -

   for( int i = Nm + 1; i < count ; i++ )
    if( tempI[ i ] > elem )
     tempI[ i ]--;

   // update the position of rho column- - - - - - - - - - - - - - - - - - -

   if( comp_col[ 0 ] > elem )
    comp_col[ 0 ]--;

   // update the items' dictionary- - - - - - - - - - -  - - - - - - - - - -

   for( Index i = 0 ; i < item_maxname ; i++ )
    if( ( dict_item[ i ] != InINF ) && ( dict_item[ i ] > elem ) )
     dict_item[ i ]--;

   // update the slack and stabilization dictionaries- - - - - - - - - - - -

   for( Index i = 0 ; i < CrrSGLen ; i++ ) {
    if( ( dict_slack[ i ] != InINF ) && ( dict_slack[ i ] > elem ) )
     dict_slack[ i ]--;
    if( stab == quadratic )
     if( ( dict_stab[ i ] != InINF ) && ( dict_stab[ i ] > elem ) )
      dict_stab[ i ]--;
    }
   }  // end for- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // compact the vocabulary - - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Compact( Upper , whch , CrrSGLen );
  Compact( Lower , whch , CrrSGLen );
  Compact( dict_slack , whch , CrrSGLen );
  if( dict_stab )
   Compact( dict_stab , whch , CrrSGLen );

  CrrSGLen = newCrrSGLen;
  }
 else {  // removing all variables- - - - - - - - - - - - - - - - - - - - - -
         // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // remove all rows in the dynamic part- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index i = 0 ; i < CrrSGLen ; i++ )
   tempI[ i ] = comp_row[ NrFi ] + i;

  osiSlvr->deleteRows( CrrSGLen , tempI );

  // mark the columns to delete - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( dict_stab )
   resizeI( 3 * CrrSGLen + item_maxname );
  else
   resizeI( 2 * CrrSGLen + item_maxname );

  count = 0;
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   if( dict_slack[ i ] != InINF ) {
    tempI[ count++ ] = dict_slack[ i ];
    if( ( Lower[ i ] > -LMINF ) &&  ( Upper[ i ] < LMINF ) )
     tempI[ count++ ] = dict_slack[ i ] + 1;
    }
   if( dict_stab )
    tempI[ count++ ] = dict_stab[ i ];
   }

  for( Index Nm = 0 ; Nm < item_maxname ; Nm++ )
   if( dict_item[ Nm ] != InINF )
    tempI[ count++ ] = dict_item[ Nm ];

  // delete the columns of all the items- - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  osiSlvr->deleteCols( count , tempI );

  // update the 0-th component- - - - - - - - - - - - - - - - - - - - - - - -

  RhoColBDm = RhoColSgPos;

  for( Index Nm = 0 ; Nm < Index( count ) ; ++Nm )
   if( comp_col[ 0 ] > tempI[ count ] )
    comp_col[ 0 ]--;

  } // end removing all elements - - - - - - - - - - - - - - - - - - - - - -
 }  // end( OSIMPSolver::RmvVars )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgAlfa( cHpRow DeltaAlfa )
{
 // note: the obj coefficients of each item column are *the opposite* of
 //       the corresponding linearization error / RHS

 MSG( 0 , "OSIMPSolver::ChgAlfa( DeltaAlfa ) called\n" );

 int numc = osiSlvr->getNumCols();
 resizeHP( numc );
 VectAssign( tempHP , osiSlvr->getObjCoefficients() , numc );

 for( Index i = 0 ; i < item_maxname ; i++ )
  if( IsSubG( i ) )
   tempHP[ dict_item[ i ] ] -= DeltaAlfa[ WComp( i ) ];
   // note the sign: the linearization error has to be *increased* by
   // DeltaAlfa[ WComp( i ) ], but since the objective function coefficients
   // are the *opposite* of the linearization error, one has to subtract

 osiSlvr->setObjective( tempHP );

 }  // end( OSIMPSolver::ChgAlfa( cHpRow ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgAlfa( cHpRow NewAlfa , cIndex wFi )
{
 // note: the obj coefficients of each item column are *the opposite* of
 //       the corresponding linearization error / RHS

 MSG( 0 , "OSIMPSolver::ChgAlfa( NewAlfa ) called\n" );

 if( wFi == 0 )
  throw( NDOException( "OSIMPSolver::ChgAlfa( * , 0 ) called" ) );

 if( wFi > NrFi ) {
  for( Index i = 0 ; i < item_maxname ; i++ )
   if( dict_item[ i ] < InINF )
    osiSlvr->setObjCoeff( dict_item[ i ] , - NewAlfa[ i ] );
  }
 else
  for( Index i = 0 ; i < item_maxname ; i++ )
   if( WComponent( i ) == wFi )
    osiSlvr->setObjCoeff( dict_item[ i ] , - NewAlfa[ i ] );

 }  // end( OSIMPSolver::ChgAlfa( cHpRow , cIndex ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgAlfa( cIndex i , cHpNum Ai )
{
 // note: the obj coefficients of each item column are *the opposite* of
 //       the corresponding linearization error / RHS

 MSG( 0 , "OSIMPSolver::ChgAlfa( i ) called\n" );

 if( ( i >= item_maxname ) || ( dict_item[ i ] == InINF ) )
  throw( NDOException( "OSIMPSolver::ChgAlfa( i ): invalid index" ) );
 
 osiSlvr->setObjCoeff( dict_item[ i ] , - Ai );

 }  // end( OSIMPSolver::ChgAlfa( cIndex , cHpNum ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChangeCurrPoint( cLMRow DLambda , cHpRow DFi )
{
 /* Note: the formula for updating the linearization error of subgradient
          i belonging to component k is

     newalfa[ i ] = alfa[ i ] + DFi[ k ] - DLambda * z[ i ] .

    (and that of a constraint is the same with DFi[ k ] replaced with 0).
    However, in OSIMPSolver:

    - the linearization error / RHS of the column corresponding to z[ i ]
      is *the opposite* of the corresponding objective function coefficient

    - the column itself contains - z[ i ]

    so beware of the sign. */

 MSG( 0 , "OSIMPSolver::ChangeCurrPoint( DLabda ) called\n" );

 // change the coefficients of the objective function- - - - - - - - - - - -

 cIndex numcols = osiSlvr->getNumCols();
 resizeHP( numcols );
 VectAssign( tempHP , osiSlvr->getObjCoefficients() , numcols );

 // update the linearization error/rhs of the items- - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < item_maxname ; i++ )
  if( dict_item[ i ] < InINF ) {
   double temp = - tempHP[ dict_item[ i ] ];
   CoinShallowPackedVector item =
                 ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );
   int k;
   for( Index j = item.getNumElements() ; j-- > 0 ; ) {
    if( ( k = item.getIndices()[ j ] - int( comp_row[ NrFi ] ) ) >= 0 )
     temp += item.getElements()[ j ] * DLambda[ k ];
    }

   if( IsSubG( i ) )
    temp += DFi[ WComp( i ) ];
   else {
    // it is a constraint: ensure its RHS remains non-negative, as small
    // negative RHS which may crop up by numerical errors would make
    // the primal master problem unfeasible (the dual unbounded)
    if( temp < 0 )
     temp = 0;
    }

   tempHP[ dict_item[ i ] ] = - temp;
   }

 // update the global lower bound (if any) - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( LwrBnds[ 0 ] > - HpINF ) {
  tempHP[ comp_col[ 0 ] ] += DFi[ 0 ];
  LwrBnds[ 0 ] -= DFi[ 0 ];
  }

 // change the cost of the easy components and the lower bound of the- - - -
 // difficult ones - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the costs of the easy component i are
 //
 //     c'[ h ] = c[ h ] - Lambda * A[ h ]
 //
 // when Lambda is updated as Lambda = Lambda + DLambda, clearly
 //
 //     c'[ h ] = c[ h ] - ( Lambda + DLambda ) * A[ h ]
 //             = ( c[ h ] - Lambda * A[ h ] ) - DLambda * A[ h ]
 //             = c'[ h ] - DLambda * A[ h ]
 
 for( Index i = 0 ; ++i <= NrFi ; )
  if( weasy[ i ] ) {  // and easy component
   cIndex BNC = FIO->GetBNC( i );
   for( Index j = 0 ; j < BNC ; j++ ) {
    double temp = tempHP[ comp_col[ i ] + j ];
    CoinShallowPackedVector item =
             ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ i ] + j );

    for( Index l = item.getNumElements() ; l > 0 ; l-- ) {
     int k;
     if( ( k = item.getIndices()[ l - 1 ] - int( comp_row[ NrFi ] ) ) >= 0 )
      temp -= item.getElements()[ l - 1 ] * DLambda[ k ];
     }

    tempHP[ comp_col[ i ] + j ] = temp;
    }
   }
  else                // a difficult component
   if( LwrBnds[ i ] > - HpINF ) {  // ... with a finite lower bound
    // the lower bound is the opposite of the linearization error and equal
    // to the coefficient: since DFi[ i ] is added to the linearization
    // error it is subtracted from the lower bound = coefficient
    tempHP[ comp_col[ i ] ] -= DFi[ i ];
    LwrBnds[ i ] -= DFi[ i ];
    }
   else                            // the lower bound is -INF
    if( ( ! NSubG[ i ] ) && ( LwrBnds[ 0 ] > - HpINF ) )  {
     // ... but the component is empty and a global lower bound is set
     // then, a "fake" individual lower bound is set to make the
     // component bounded, whose value has to be < the global lower bound
     tempHP[ comp_col[ i ] ] = LwrBnds[ 0 ] - 1;
     }

 // change the bounds of the variables - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < CrrSGLen ; i++ ) {
  if( Upper[ i ] < LMINF )
   Upper[ i ] -= DLambda[ i ];
  if( Lower[ i ] > -LMINF )
   Lower[ i ] -= DLambda[ i ];
  }

 // update the slack prices - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ptUpdatePricesInPlace();

 // finally, change all costs in one blow - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 osiSlvr->setObjective( tempHP );

 }  // end( OSIMPSolver::ChangeCurrPoint( DLambda , DFi ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChangeCurrPoint( cHpNum Tau , cHpRow DFi )
{
 MSG( 0 , "OSIMPSolver::ChangeCurrPoint( Tau ) called\n" );

 cIndex numcols = osiSlvr->getNumCols();
 resizeHP( numcols );
 VectAssign( tempHP , osiSlvr->getObjCoefficients() , numcols );

 #if( PRESERVE_OSI_SOLS == 0 )
  const double *rsol = osiSlvr->getRowPrice();
 #endif
 const double *delta = rsol + comp_row[ NrFi ];

 // update the linearization error/rhs of the items- - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < item_maxname ; i++ )
  if( dict_item[ i ] < InINF ) {
   double temp = - tempHP[ dict_item[ i ] ];

   temp -= ( Tau / t ) * GiPerd[ i ];

   if( IsSubG( i ) )
    temp += DFi[ WComp( i ) ];
   else {  // it is a constraint: ensure its RHS remains non-negative, as
    // small negative RHS which may crop up by numerical errors may make
    // the primal master problem unfeasible (the dual unbounded)
    if( temp < 0 )
     temp = 0;
    }

   tempHP[ dict_item[ i ] ] = - temp;
   }

 // update the global lower bound (if any) - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( LwrBnds[ 0 ] > - HpINF ) {
  tempHP[ comp_col[ 0 ] ] += DFi[ 0 ];
  LwrBnds[ 0 ] -= DFi[ 0 ];
  }

 // change the cost of the easy components and the lower bound of the- - - -
 // difficult ones - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; ++i <= NrFi ; )
  if( weasy[ i ] ) {  // an easy component
   if( Tau > 0 ) {  // when Tau == 0, the costs do not change
    cIndex BNC = FIO->GetBNC( i );
    for( Index j = 0 ; j < BNC ; j++ ) {
     double temp = tempHP[ comp_col[ i ] + j ];
     CoinShallowPackedVector item =
             ( osiSlvr->getMatrixByCol() )->getVector( comp_col[ i ] + j );

     for( Index l = item.getNumElements() ; l > 0 ; l-- ) {
      int k;
      if( ( k = item.getIndices()[ l - 1 ] - int( comp_row[ NrFi ] ) ) >= 0 )
       temp -= item.getElements()[ l - 1 ] * delta[ k ] * ( Tau / t );
      }

     tempHP[ comp_col[ i ] + j ] = temp;
     }
    }
   }
  else                // a difficult component
   if( LwrBnds[ i ] > - HpINF ) {  // ... with a finite lower bound
    // the lower bound is the opposite of the linearization error and equal
    // to the coefficient: since DFi[ i ] is added to the linearization
    // error it is subtracted from the lower bound = coefficient
    tempHP[ comp_col[ i ] ] -= DFi[ i ];
    LwrBnds[ i ] -= DFi[ i ];
    }
   else                            // the lower bound is -INF
    if( ( ! NSubG[ i ] ) && ( LwrBnds[ 0 ] > - HpINF ) )  {
     // ... but the component is empty and a global lower bound is set
     // then, a "fake" individual lower bound is set to make the
     // component bounded, whose value has to be < the global lower bound
     tempHP[ comp_col[ i ] ] = LwrBnds[ 0 ] - 1;
     }

 // change the bounds of the variables and therefore the slack prices- - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( Tau > 0 ) {  // when Tau == 0, nothing changes
  cLMRow d = Readd( true );
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   if( Upper[ i ] < LMINF )
    Upper[ i ] -= d[ i ] * ( Tau / t );
   if( Lower[ i ] > -LMINF )
    Lower[ i ] -= d[ i ] * ( Tau / t );
   }

  ptUpdatePricesInPlace();
  }

 // finally, change all costs in one blow - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 osiSlvr->setObjective( tempHP );

 }  // end( OSIMPSolver::ChangeCurrPoint( Tau , DFi ) )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgSubG( cIndex strt , Index stp , cIndex wFi )
{
 MSG( 0 , "OSIMPSolver::ChgSubG() called\n" );

 Index str_p = comp_row[ NrFi ];
 Index count = osiSlvr->getNumCols();

 std::vector< double > lower;
 std::vector< double > upper;

 if( NrEasy > 1 ) {
  lower.resize( count );
  VectAssign( lower.data() , osiSlvr->getColLower() , count );
  upper.resize( count );
  VectAssign( upper.data() , osiSlvr->getColUpper() , count );
  }
 
 resizeHP( str_p + CrrSGLen );
 resizeI( str_p + CrrSGLen );

 std::vector< int > ColsToDelete;

 // the 0-th component of Fi has changed - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ! wFi ) || ( wFi == InINF ) ) {
  // get the item  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // to make life simpler get it from strt to the end and rebuild it all

  cIndex_Set SGBse;
  Index SGBDim = FIO->GetGi( tempHP , SGBse , FIO->GetMaxName() ,
			     strt , CrrSGLen );

  // find the position in the "subgradient part" of the sparse RhoCol[] of
  // the first index >= strt; note that the first row of the "subgradient
  // part" is comp_row[ NrFi ], and therefore the index value to look for
  // is comp_row[ NrFi ] + strt
  Index i = RhoColSgPos;  // only change the 
  while( ( i < RhoColBDm ) && ( RhoColBse[ i ] < comp_row[ NrFi ] + strt ) )
   i++;

  // once the position is found, reset the "subgradient part" of the sparse
  // RhoCol[] from there on and rebuild it anew
  RhoColBDm = i;
  if( SGBse ) {
   for( Index j = 0 ; j < SGBDim ; ++j )
    if( tempHP[ j ] ) {
     RhoColBse[ RhoColBDm ] = comp_row[ NrFi ] + SGBse[ j ];
     RhoCol[ RhoColBDm++ ] = - tempHP[ j ];
     }
   }
  else
   for( Index j = 0 ; j < SGBDim ; ++j )
    if( tempHP[ j ] ) {
     RhoColBse[ RhoColBDm ] = comp_row[ NrFi ] + j;
     RhoCol[ RhoColBDm++ ] = - tempHP[ j ];
     }

  UpdateRhoCol();  // remove the previous version and add the new one
  }

 // only the wFi-th component of Fi has changed- - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( 1 <= wFi ) && ( wFi <= NrFi ) ) {

  if( weasy[ wFi ] ) {  // if the component is an easy one - - - - - - - - - -
                        // - - - - - - - - - - - - - - - - - - - - - - - - - -

   // allocate the memory to describe the matrix A[ i ]- - - - - - - - - - - -

   int *Abeg = new int[ FIO->GetBNC( wFi ) + 1 ];
   int *Aind = new int[ FIO->GetANZ( wFi ) ];
   double *Aval = new double[ FIO->GetANZ(  wFi ) ];

   int col = comp_col[ wFi ];  // get the name of the column
   comp_col[ wFi ] = osiSlvr->getNumCols();
   for( Index i = 0 ; i < FIO->GetBNC( wFi ) ; i++ ) {

    // copy the unchanged part of the item, if exists- - - - - - - - - - - - -

    CoinShallowPackedVector p_item =
        	     ( osiSlvr->getMatrixByCol() )->getVector( col + i );

    count = 0;
    for( int j = 0 ; j < p_item.getNumElements() ; j++ )
     if( p_item.getIndices()[ j ] < ( int( str_p ) + strt ) ) {
      tempI[ count ] = p_item.getIndices()[ j ];
      tempHP[ count++ ] = p_item.getElements()[ j ];
      }
     else
      if( p_item.getIndices()[ j ] >= int( str_p ) + stp ) {
       tempI[ count ] = p_item.getIndices()[ j ];
       tempHP[ count++ ] = p_item.getElements()[ j ];
       }

    // ask FiOracle for the matrix A[ i ], it may give a partial
    // information of the rows of matrix A[ i ]- - - - - - - - - - - - - - - -

    FIO->GetADesc( wFi , Abeg , Aind , Aval , strt , stp );

    // write the dynamic part: A's columns - - - - - - - - - - - - - - - - - -

    for( int k = Abeg[ i ] ; k < Abeg[ i + 1 ] ; k++ ) {
     tempI[ count ] = Aind[ k ] + str_p;
     tempHP[ count++ ] = Aval[ k ];
     }

    //  ... add the cost
    HpNum coeff = osiSlvr->getObjCoefficients()[ col + i ];

    // set the item- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ColsToDelete.push_back( col + i );
    osiSlvr->addCol( count , tempI , tempHP , lower[ col + i ] ,
		     upper[ col + i ] , coeff );
    }

   delete[] Abeg;
   delete[] Aind;
   delete[] Aval;
   }
  else // if wFi is a difficult component
   for( Index i = 0 ; i < item_maxname ; i++ )
    if( WComponent( i ) == wFi ) {
     // copy the unchanged part of the item, if exists - - - - - - - - - - - -

     CoinShallowPackedVector p_item =
    	     ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

     count = 0;
     for( int j = 0 ; j < p_item.getNumElements() ; j++ )
      if( p_item.getIndices()[ j ] < ( int( str_p ) + strt ) ) {
       tempI[ count ] = p_item.getIndices()[ j ];
       tempHP[ count++ ] = p_item.getElements()[ j ];
       }
      else
       if( p_item.getIndices()[ j ] >= int( str_p ) + stp ) {
        tempI[ count ] = p_item.getIndices()[ j ];
        tempHP[ count++ ] = p_item.getElements()[ j ];
        }

     // get the item - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     cIndex_Set SGBse;
     Index SGBDim = FIO->GetGi( tempHP + count , SGBse , i , strt , stp );
     VectScale( tempHP + count , - double( 1 ) , SGBDim  );

     if( SGBse == nullptr )
      std::iota( tempI + count , tempI + count + SGBDim , str_p + strt );
     else
      for( Index j = 0 ; j < SGBDim ; j++ )
       tempI[ j + count ] =  SGBse[ j ] + str_p;

     SGBDim += count;

     int col = dict_item[ i ];  // get the name of the column
     HpNum coeff = osiSlvr->getObjCoefficients()[ col ]; // .. and its price

     // set the item - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     ColsToDelete.push_back( col );
     dict_item[ i ] = osiSlvr->getNumCols();
     osiSlvr->addCol( SGBDim , tempI , tempHP , 0 , osiSlvr->getInfinity() ,
		      coeff );
     }
  }
 else
  if( wFi > NrFi  ) {  // all the components of Fi have changed

   for( Index wFiComp = 1 ; wFiComp <= NrFi ; wFiComp++ )
    if( weasy[ wFiComp ] ) {  // if the component is an easy one

     // allocate the memory to describe the matrix A[ i ]
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     int *Abeg = new int[ FIO->GetBNC( wFiComp ) + 1 ];
     int *Aind = new int[ FIO->GetANZ( wFiComp ) ];
     double *Aval = new double[ FIO->GetANZ( wFiComp ) ];

     int col = comp_col[ wFiComp ];  // get the name of the column
     comp_col[ wFiComp ] = osiSlvr->getNumCols();
     for( Index i = 0 ; i < FIO->GetBNC( wFiComp ) ; i++ ) {

      // copy the unchanged part of the item, if exists  - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      CoinShallowPackedVector p_item =
        	     ( osiSlvr->getMatrixByCol() )->getVector( col + i );

      count = 0;
      for( int j = 0 ; j < p_item.getNumElements() ; j++ )
       if( p_item.getIndices()[ j ] < ( int( str_p ) + strt ) ) {
        tempI[ count ] = p_item.getIndices()[ j ];
        tempHP[ count++ ] = p_item.getElements()[ j ];
        }
       else
        if( p_item.getIndices()[ j ] >= int( str_p ) + stp ) {
         tempI[ count ] = p_item.getIndices()[ j ];
         tempHP[ count++ ] = p_item.getElements()[ j ];
         }

      // ask FiOracle for the matrix A[ i ], it  may give a partial
      // information of the rows of matrix A[ i ]
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      FIO->GetADesc( wFiComp , Abeg , Aind , Aval , strt , stp );

      // write the dynamic part: A's columns - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      for( int k = Abeg[ i ] ; k < Abeg[ i + 1 ] ; k++ ) {
       tempI[ count ] = Aind[ k ] + str_p;
       tempHP[ count++ ] = Aval[ k ];
       }

      // .. add the cost
      HpNum coeff = osiSlvr->getObjCoefficients()[ col + i ];

      // set the item  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


      ColsToDelete.push_back( col + i );
      osiSlvr->addCol( count , tempI , tempHP , lower[ col + i ] ,
		       upper[ col + i ] , coeff );
      }

     delete[] Abeg;
     delete[] Aind;
     delete[] Aval;
     }

   for( Index i = 0 ; i < item_maxname ; i++ )
    if( WComponent( i ) < InINF ) {

     // copy the unchanged part of the item, if exists   - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     CoinShallowPackedVector p_item =
      ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

     count = 0;
     for( int j = 0 ; j < p_item.getNumElements() ; j++ )
	  if( p_item.getIndices()[ j ] < ( int( str_p ) + strt ) ) {
       tempI[ count ] = p_item.getIndices()[ j ];
       tempHP[ count++ ] = p_item.getElements()[ j ];
       }
      else
       if( p_item.getIndices()[ j ] >= int( str_p ) + stp ) {
        tempI[ count ] = p_item.getIndices()[ j ];
        tempHP[ count++ ] = p_item.getElements()[ j ];
        }

     // get the item   - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     cIndex_Set SGBse;
     Index SGBDim = FIO->GetGi( tempHP + count , SGBse , i , strt , stp );
     VectScale( tempHP + count , -double(1) , SGBDim  );

     if( SGBse == nullptr )
      std::iota( tempI + count , tempI + count + SGBDim , str_p + strt );
     else
      for( Index j = 0 ; j < SGBDim ; j++ )
	   tempI[ j + count ] =  SGBse[ j ] + str_p;

     SGBDim += count;

     int col = dict_item[ i ];  // get the name of the column
     HpNum coeff = osiSlvr->getObjCoefficients()[ col ]; // .. and add price

     // set the item    - - - - - - - - - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     ColsToDelete.push_back( col );
     dict_item[ i ] = osiSlvr->getNumCols();
     osiSlvr->addCol( SGBDim , tempI , tempHP , 0 , osiSlvr->getInfinity() ,
		      coeff );
     }
   }

 // update the dictionaries   - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::sort( ColsToDelete.begin() , ColsToDelete.end() , std::greater< int >() );
 for( auto & el : ColsToDelete ) {
  for( Index i = 0 ; i < item_maxname ; i++ )
   if( ( dict_item[ i ] < InINF ) && ( dict_item[ i ] > el ) )
    dict_item[ i ]--;

  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   if( ( dict_slack[ i ] < InINF ) && ( dict_slack[ i ] > el ) )
    dict_slack[ i ]--;

   if( dict_stab && ( dict_stab[ i ] < InINF ) && ( dict_stab[ i ] > el ) )
    dict_stab[ i ]--;
   }

  for( Index i = 0 ; i <= NrFi ; i++ )
   if( ( comp_col[ i ] > el ) && ( comp_col[ i ] < InINF ) )
    comp_col[ i ]--;
  }

 osiSlvr->deleteCols( ColsToDelete.size() , ColsToDelete.data() );

 }  // end( OSIMPSolver::ChgSubG )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgCosts( Index wFi , cLMRow Lambda )
{
 if( ( ! weasy ) || ( ! weasy[ wFi ] ) )
  throw( std::invalid_argument( "ChgCosts called on non-easy component" ) );

 auto BNC = FIO->GetBNC( wFi );
 /*!! here one should check that the number of columns has not changed, but
      since variables can be moved around there is no easy way of knowing
      how many columns were there initially
 if( BNC != comp_col[ wFi + 1 ] - comp_col[ wFi ] )
  throw( std::invalid_argument(
		       "ChgCosts: changing number of columns not allowed" ) );
		       !!*/

 // allocate memory for the costs
 auto cst = new double[ BNC ];

 // get the costs from the FiOracle
 FIO->GetBDesc( wFi , nullptr , nullptr , nullptr , nullptr , nullptr , cst ,
		nullptr , nullptr );

 // the costs to be set in the master problem are the Lagrangian ones
 // c[ i ] - Lambda * A[ i ], hence retrieve A[ i ]

 auto ANZ = FIO->GetANZ( wFi );
 auto Abeg = new int[ BNC + 1 ];
 auto Aind = new int[ ANZ ];
 auto Aval = new double[ ANZ ];

 FIO->GetADesc( wFi , Abeg , Aind , Aval );

 // compute and set the Lagrangian costs
 for( Index j = 0 ; j < BNC ; ++j ) {
  auto cj = cst[ j ];
   
  for( int k = Abeg[ j ] ; k < Abeg[ j + 1 ] ; ++k )
   cj -= Aval[ k ] * Lambda[ Aind[ k ] ];

  osiSlvr->setObjCoeff( comp_col[ wFi ] + j , cj );
  }

 delete[] Aval;
 delete[] Aind;
 delete[] Abeg;

 delete[] cst;

 }  // end( OSIMPSolver::ChgCosts )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgRLHS( Index wFi )
{
 if( ( ! weasy ) || ( ! weasy[ wFi ] ) )
  throw( std::invalid_argument( "ChgRLHS called on non-easy component" ) );

 auto BNR = FIO->GetBNR( wFi );
 if( ! BNR )  // nothing to change
  return;     // all done

 auto BNC = FIO->GetBNC( wFi );
 /*!! here one should check that the number of columns has not changed, but
      since variables can be moved around there is no easy way of knowing
      how many columns were there initially
 if( BNC != comp_col[ wFi + 1 ] - comp_col[ wFi ] )
  throw( std::invalid_argument(
		       "ChgRLHS: changing number of columns not allowed" ) );
		       ||*/

 auto lbd = new double[ BNC ];
 auto ubd = new double[ BNC ];
 auto lhs = new double[ BNR ];
 auto rhs = new double[ BNR ];

 FIO->GetBDesc( wFi , 0 , 0 , 0 , lhs , rhs , 0 , lbd , ubd );

 // find the starting position of the rows of wFi
 const auto stp = RhoColBse + RhoColSgPos;
 auto it = std::lower_bound( RhoColBse , stp , comp_row[ wFi - 1 ] );
 auto rit = RhoCol + std::distance( RhoColBse , it );
 
 // first skip the (unchanged) part about upper/lower bounds
 for( Index j = 0 ; j < BNC ; ++j ) {
  if( lbd[ j ] && ( lbd[ j ] > -Inf< double >() ) )
   // l[ i ] <= x[ i ] becomes  0 <= x[ i ] - \rho l[ i ] (< INF)
   rit++;

  if( ubd[ j ] && ( ubd[ j ] < Inf< double >() ) )
   // x[ i ] <= u[ i ] becomes  (-INF) <= x[ i ] - \rho u[ i ] <= 0
   rit++;
  }

 // now change the part about LHS/RHS: note that we assume that row
 // numbers remain the same, as they should
 for( Index j = 0 ; j < BNR ; ++j ) {  // for each row j
  // if both lhs and rhs are non-0 and non-inf, two rows are created
  if( lhs[ j ] && ( lhs[ j ] > -Inf< double >() ) &&
      rhs[ j ] && ( rhs[ j ] < Inf< double >() ) &&
      ( lhs[ j ] < rhs[ j ] ) ) {
   *(rit++) = - lhs[ j ];  // hence RhoCol has a nonzero for the lhs
   *(rit++) = - rhs[ j ];  // ... and a separate nonzero for the rhs
   }
  else {                           // one row is created
   double val = 0;                 // ... but it may have a 0 RHS
   if( lhs[ j ] && ( lhs[ j ] > -Inf< double >() ) )
    val = - lhs[ j ];
   else
    if( rhs[ j ] && ( rhs[ j ] < Inf< double >() ) )
     val = - rhs[ j ];

   // a nonzero in RhoCol is only created if one (and only one) among the
   // rhs and lhs is finite and nonzero; note that one single row still have
   // been created for this constraint, but there is no nonzero coefficient
   // for that row in RhoCol (the assumption being exactly that since it's
   // not here now it has never been there before)
   if( val )
    *(rit++) = val;
   }
  }

 delete[] rhs;
 delete[] lhs;
 delete[] ubd;
 delete[] lbd;

 /* This check had originally been written with "!=" instead of ">", but in
  * that form it was wrong in that is assumed that the number of nonzeros in
  * RhoCol for this component is the same as the number of rows in the master
  * problem for this component. This is *not* true in general because rows in
  * the master problem that naturally have only 0 or +/-INF lhs/rhs do not
  * produce a nonzero in RhoCol; hence, the number of nonzero in RhoCol must
  * be <= than the number of rows, but it can be < (and therefore !=).
  * Note that this only concerns "real" rows, as bound constraints that have
  * only 0 or +/-INF lbd/ubd do not produce a row in the master problem
  * (because they rather produce a bound constraint). */

 if( std::distance( RhoCol , rit ) - std::distance( RhoColBse , it ) >
     comp_row[ wFi ] - comp_row[ wFi - 1 ] )
  throw( std::invalid_argument(
		          "ChgRLHS: changing number of rows not allowed" ) );

 
 UpdateRhoCol();  // remove the previous version and add the new one

 }  // end( OSIMPSolver::ChgRLHS )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ChgLUBD( Index wFi )
{
 if( ( ! weasy ) || ( ! weasy[ wFi ] ) )
  throw( std::invalid_argument( "ChgLUBD called on non-easy component" ) );

 auto BNC = FIO->GetBNC( wFi );
 /*!! here one should check that the number of columns has not changed, but
      since variables can be moved around there is no easy way of knowing
      how many columns were there initially
 if( BNC != comp_col[ wFi + 1 ] - comp_col[ wFi ] )
  throw( std::invalid_argument(
		        "ChgLUBD: changing number of columns not allowed" ) );
			!!*/

 auto lbd = new double[ BNC ];
 auto ubd = new double[ BNC ];

 FIO->GetBDesc( wFi , nullptr , nullptr , nullptr , nullptr , nullptr ,
		nullptr , lbd , ubd );

 const auto stp = RhoColBse + RhoColSgPos;
 auto it = std::lower_bound( RhoColBse , stp , comp_row[ wFi - 1 ] );
 auto rit = RhoCol + std::distance( RhoColBse , it );

 Index boundrow = comp_row[ wFi - 1 ];   // first row of bounds replacement

 for( Index j = 0 ; j < BNC ; ++j ) {
  double lbj = lbd[ j ];
  if( lbj && ( lbj > -Inf< double >() ) ) {
   // l[ i ] <= x[ i ] becomes  0 <= x[ i ] - \rho l[ i ] (< INF)
   if( *(it++) != boundrow++ )
    throw( std::logic_error( "change in bound structure not allowed" ) );
   *(rit++) = - lbj;
   lbj = - osiSlvr->getInfinity();
   }

  double ubj = ubd[ j ];
  if( ubj && ( ubj < Inf< double >() ) ) {
   // x[ i ] <= u[ i ] becomes  (-INF) <= x[ i ] - \rho u[ i ] <= 0
   if( *(it++) != boundrow++ )
    throw( std::logic_error( "change in bound structure not allowed" ) );
   *(rit++) = - ubj;
   ubj = osiSlvr->getInfinity();
   }

  osiSlvr->setColBounds( comp_col[ wFi ] + j , lbj , ubj );
  }

 delete[] ubd;
 delete[] lbd;

 UpdateRhoCol();  // remove the previous version and add the new one

 }  // end( OSIMPSolver::ChgLUBD )

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

OSIMPSolver::~OSIMPSolver()
{
 MSG( 0 , "OSIMPSolver::~OSIMPSolver() called\n" );

 cleanup();

 delete derhand;

 }  // end( OSIMPSolver::~OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void OSIMPSolver::cleanup( void )
{
 if( osiSlvr )
  osiSlvr->reset();

 // delete the dictionaries everything

 delete[] tempHP;
 tempHP = nullptr;
 tempHP_size = 0;
 delete[] tempI;
 tempI = nullptr;
 tempI_size = 0;

 delete[] dict_stab;
 dict_stab = nullptr;

 delete[] GiPerd;
 GiPerd = nullptr;

 delete[] LwrBnds;
 LwrBnds = nullptr;

 delete[] RhoCol;
 RhoCol = nullptr;
 delete[] RhoColBse;
 RhoColBse = nullptr;
 RhoColBDm = RhoColSgPos = 0;

 delete[] Lower;
 Lower = nullptr;
 delete[] Upper;
 Upper = nullptr;

 delete[] wcomp;
 wcomp = nullptr;

 delete[] dict_slack;
 dict_slack = nullptr;
 delete[] dict_item;
 dict_item = nullptr;

 delete[] comp_col;
 comp_col = nullptr;
 delete[] comp_row;
 comp_row = nullptr;
 item_maxname = 0;

 delete[] NConst;
 NConst = nullptr;
 delete[] NSubG;
 NSubG = nullptr;

 delete[] NewItem;
 NewItem = nullptr;

 delete[] weasy;
 weasy = nullptr;

 #if( PRESERVE_OSI_SOLS )
  delete[] rcst;
  rcst = nullptr;
  delete[] rsol;
  rsol = nullptr;
  rsols = 0;
  delete[] csol;
  csol = nullptr;
  csols = 0;
 #endif

 t = 1;
 first = true;
 MaxBSize = 0;
    
 NRCall = 0;          

 } // end( OSIMPSolver::cleanup )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::tUpdatePrices( cIndex strt , Index stp )
{
 if( ! dict_slack )  // no slacks
  return;            // nothing to do

 if( stp > CrrSGLen )
  stp = CrrSGLen;

 if( stab == boxstep ) {
   resizeI( 2 * ( stp - strt ) );
   resizeHP( 2 * ( stp - strt ) );
   int *ttI = tempI;
   HpRow ttHP = tempHP;

   for( Index i = strt ; i < stp ; i++ ) {
    *(ttI++) = dict_slack[ i ];
    *(ttHP++) = ( Lower[ i ] > -t ) ? Lower[ i ] : -t;
    *(ttI++) = dict_slack[ i ] + 1;
    *(ttHP++) = ( t > Upper[ i ] ) ? -Upper[ i ] : -t;
    }

   osiSlvr->setObjCoeffSet( tempI , ttI , tempHP );
   }
 else
  if( stab == quadratic ) {

   #if WHICH_OSI_MP == 1
    resizeHP( osiSlvr->getNumCols() );
    VectAssign( tempHP , HpNum( 0 ) , osiSlvr->getNumCols() );
    for( Index i = 0 ; i < CrrSGLen ; ++i )
     tempHP[ dict_stab[ i ] ] = -t;

    auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
    CPXcopyqpsep( osiCpx->getEnvironmentPtr() , osiCpx->getLpPtr() , tempHP );
   #elif WHICH_OSI_MP == 2
    // note that the primal stabilising term is ( t / 2 ) || s ||_2^2
    // Gurobi does *not* automatically add the "standard" 1 / 2
    // coefficient in front of the quadratic term, while other solvers
    // (namely, Cplex) do, which requires to explicitly scale t here
    resizeHP( CrrSGLen );
    VectAssign( tempHP , -t / 2 , CrrSGLen );
    resizeI( CrrSGLen );
    for( Index i = 0 ; i < CrrSGLen ; ++i )
     tempI[ i ] = dict_stab[ i ];

    auto osiGrb = dynamic_cast< OsiGrbSolverInterface * >( osiSlvr );
    auto GrbLp = osiGrb->getLpPtr();
    // note that GRBaddqpterms() *adds* the new quadratic terms to the
    // existing ones, which is why one has to call GRBdelq() before in
    // order to zero the quadratic part of the objective and renew it
    GRBdelq( GrbLp );
    GRBaddqpterms( GrbLp , CrrSGLen , tempI , tempI , tempHP );
   #else
    throw( NDOException( "OSIMPSolver: stab == quadratic not supported" ) );
   #endif
   }
 }  // end( OSIMPSolver::tUpdatePrices )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ptUpdatePrices( cIndex strt , Index stp )
{
 // distinguish box stabilization case because it depends on t too- - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! dict_slack )  // no slacks
  return;            // nothing to do

 if( stp > CrrSGLen )
  stp = CrrSGLen;

 resizeI( 2 * ( stp - strt ) );
 resizeHP( 2 * ( stp - strt ) );
 int *ttI = tempI;
 HpRow ttHP = tempHP;

 if( stab != boxstep )  // anything but boxstep
  for( Index i = strt ; i < stp ; i++ ) {
   int offset = 0;  // if any lower bound exists, there is a slack s_i^+
   if( Lower[ i ] > -LMINF ) {
    *(ttI++) = dict_slack[ i ];
    *(ttHP++) = Lower[ i ];
    offset = 1;
    }
   if( Upper[ i ] < LMINF ) {
    *(ttI++) = dict_slack[ i ] + offset;
    *(ttHP++) = -Upper[ i ];
    }
   }
 else                   // boxstep
  for( Index i = strt ; i < stp ; i++ ) {
   *(ttI++) = dict_slack[ i ];
   *(ttHP++) = ( Lower[ i ] > -t ) ? Lower[ i ] : -t;
   *(ttI++) = dict_slack[ i ] + 1;
   *(ttHP++) = ( t > Upper[ i ] ) ? -Upper[ i ] : -t;
   }

 osiSlvr->setObjCoeffSet( tempI , ttI , tempHP );

 }  // end( OSIMPSolver::ptUpdatePrices )

/*--------------------------------------------------------------------------*/

void OSIMPSolver::ptUpdatePricesInPlace( void )
{
 if( ! dict_slack )  // no slacks
  return;            // nothing to do

 if( stab != boxstep )  // anything but boxstep
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   int offset = 0;  // if any lower bound exists, there is a slack s_i^+
   if( Lower[ i ] > -LMINF ) {
    tempHP[ dict_slack[ i ] ] = Lower[ i ];
    offset = 1;
    }
   if( Upper[ i ] < LMINF )
    tempHP[ dict_slack[ i ] + offset ] = -Upper[ i ];
   }
 else                   // boxstep
  for( Index i = 0 ; i < CrrSGLen ; i++ ) {
   tempHP[ dict_slack[ i ] ] = ( Lower[ i ] > -t ) ? Lower[ i ] : -t;
   tempHP[ dict_slack[ i ] + 1 ] = ( t > Upper[ i ] ) ? -Upper[ i ] : -t;
   }

 }  // end( OSIMPSolver::ptUpdatePricesInPlace )

/*------------------------------------------------------------------------*/

void OSIMPSolver::UpdateRhoCol( void )
{
 int col = comp_col[ 0 ];  // get the name of the column of rho- - - - - - -
 HpNum coeff = osiSlvr->getObjCoefficients()[ col ]; // and the lower bound

 osiSlvr->deleteCols( 1 , &col ); // delete the previous rho column

 // update the dictionaries- - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index i = 0 ; i < item_maxname ; i++ )
  if( ( dict_item[ i ] < InINF ) && ( dict_item[ i ] > Index( col ) ) )
   dict_item[ i ]--;

 for( Index i = 0 ; i < CrrSGLen ; i++ )
  if( ( dict_slack[ i ] < InINF ) && ( dict_slack[ i ] > Index( col ) ) )
   dict_slack[ i ]--;

 if( dict_stab )
  for( Index i = 0 ; i < CrrSGLen ; i++ )
   if( ( dict_stab[ i ] < InINF ) && ( dict_stab[ i ] > Index( col ) ) )
    dict_stab[ i ]--;

 for( Index i = 0 ; ++i <= NrFi ; )
  if( ( comp_col[ i ] > Index( col ) ) && ( comp_col[ i ] < InINF ) )
   comp_col[ i ]--;

 // add the rho column to the master problem - - - - - - - - - - - - - - - -

 comp_col[ 0 ] = osiSlvr->getNumCols();
 osiSlvr->addCol( RhoColBDm , RhoColBse , RhoCol , 0 ,
		  osiSlvr->getInfinity() , coeff );

 }  // end( OSIMPSolver::UpdateRhoCol )

/*------------------------------------------------------------------------*/

void OSIMPSolver::switchToQP( void )
{
 #if WHICH_OSI_MP == 1
  auto osiCpx = dynamic_cast< OsiCpxSolverInterface * >( osiSlvr );
  if( CPXchgprobtype( osiCpx->getEnvironmentPtr() ,  osiCpx->getLpPtr() ,
		      CPXPROB_QP ) )
   throw( NDOException( "OSIMPSolver::switchToQP: can't turn to QP problem"
			) );
 #elif WHICH_OSI_MP == 2
  // nothing to do for Gurobi
 #else
  throw( NDOException( "OSIMPSolver::switchToQP: not implemented yet" ) );
 #endif

 } // end( OSIMPSolver::switchToQP )

/*------------------------------------------------------------------------*/

/* bool OSIMPSolver::isactive( Index i )
{
 return( osiSlvr->getRowSense()[ comp_row[ NrFi ] + i ] == 'E' );
 } */

/*------------------------------------------------------------------------*/

void OSIMPSolver::activate( Index i )
{
 osiSlvr->setRowType( comp_row[ NrFi ] + i , 'E' , 0 , 0. );
 }

/*------------------------------------------------------------------------*/

void OSIMPSolver::deactivate( Index i )
{
 osiSlvr->setRowType( comp_row[ NrFi ] + i , 'N' , 0 , 0 );
 }

/*------------------------------------------------------------------------*/

void OSIMPSolver::resizeHP( Index i )
{
 if( i > tempHP_size ) {
  delete[] tempHP;
  tempHP = new HpNum[ tempHP_size = i ];
  }
 }

/*------------------------------------------------------------------------*/

void OSIMPSolver::resizeI( Index i )
{
 if( i > tempI_size ) {
  delete[] tempI;
  tempI = new int[ tempI_size = i ];
  }
 }

/*------------------------------------------------------------------------*/

void OSIMPSolver::CheckDS( void )
{
 #if CHECK_DS & 1
  // obvious sanity check in dict_item[]
  if( osiSlvr ) {
   Index numcols = osiSlvr->getNumCols();
   for( Index name = 0 ; name < item_maxname ; name++ )
    if( dict_item[ name ] < InINF )
     if( dict_item[ name ] > numcols )
      cout << "dict_item[ " << name << " ] = " << dict_item[ name ]
	   << " out of range (> " << numcols << ")" << std::endl;
   }
 #endif
 }

/*--------------------------------------------------------------------------*/

Index OSIMPSolver::CheckBCopy( void )
{
 // check if the item is identical to any other in the bundle - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index IsIde = InINF;

 if( checkID ) {
  if( NewItemBse ) {  // sparse format for the new item - - - - - - - - - - -
                      //- - - - - - - - - - - - - - - - - - - - - - - - - - -
   int h;
   Index NumElem;
   for( Index i = 0 ; i < item_maxname ; i++ )
    if( dict_item[ i ] < InINF  )  {
     // the items can't be of different type- - - - - - - - - - - - - - - - -
     if( NewItemisSG ^ IsSubG( i ) )
      continue;

     // the subgradients must be relative to the same component - - - - - - -
     if( NewItemisSG && ( WComp( i ) != NewItemFi ) )
      continue;

     // get the item- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     CoinPackedVector item =
      ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

     item.sortIncrIndex();
     NumElem = item.getNumElements();

     const int *B1 = item.getIndices();
     cIndex_Set B2 = NewItemBse;

     cSgRow g1 = item.getElements();
     cSgRow g2 = NewItem;

     // take off the static part of the item - - - - - - - - - - - - - - - -
     for( ; ( ( h = ( int( *B1 ) - int( comp_row[ NrFi ] ) ) ) < 0 ) &&
	    NumElem ; ) {
      B1++;
      g1++;
      NumElem--;
      }

     // checks whether or not the two items are element-wise identical - - -
     if( NumElem != NewItemBDm )
      continue;

     for( ; NumElem ; NumElem-- ) {
      if( ( *(B1++) - comp_row[ NrFi ] ) != *(B2++) )
       break;

      if( ( - *(g1++) ) != *(g2++) )
       break;
      }

     if( NumElem == 0 ) {
      IsIde = i;
      break;
      }
     }
   }
  else {  // dense format for the new item- - - - - - - - - - - - - - - - - -
	  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   int h;
   Index l;
   Index NumElem;

   for( Index i = 0 ; i < item_maxname ; i++ )
    if( dict_item[ i ] < InINF )  {
     // the items can't be of different type- - - - - - - - - - - - - - - - -
     if( NewItemisSG ^ IsSubG( i ) )
      continue;

     // the subgradients must be relative to the same component - - - - - - -
     if( NewItemisSG && ( WComp( i ) != NewItemFi ) )
      continue;

     // get the item- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     CoinPackedVector item =
      ( osiSlvr->getMatrixByCol() )->getVector( dict_item[ i ] );

     item.sortIncrIndex();
     NumElem = item.getNumElements();
     l = 0;

     const int *B = item.getIndices();

     cSgRow g1 = NewItem;
     cSgRow g2 = item.getElements();

     // take off the static part of the item- - - - - - - - - - - - - - - - -
     for( ; NumElem &&
	    ( ( h = ( int( *B ) - int( comp_row[ NrFi ] ) ) ) < 0 ) ; ) {
      B++;
      g2++;
      NumElem--;
      }

     // checks whether or not the two items are element-wise identical- - - -
     for( ; NumElem ; NumElem-- , l++ ) {
      h = *(B++) - comp_row[ NrFi ];

      for( ; l < h ; l++ )
       if( *(g1++) )
	break;

      if( l < h )
       break;

      if( *(g1++) != ( -*(g2++) ) )
       break;
      }

     // the last part of g1 should be zero- - - - - - - - - - - - - - - - - -

     if( NumElem )
      continue;

     for( ; l < CrrSGLen ; l++ )
      if( *(g1++) ) {
       NumElem = InINF;
       break;
       }

     if( NumElem == 0 ) {
      IsIde = i;
      break;
      }
     }
   }
  }

 return( IsIde );

 }  // end( OSIMPSolver )

/*--------------------------------------------------------------------------*/
/*------------------------- End File OSIMPSolver.C -------------------------*/
/*--------------------------------------------------------------------------*/
