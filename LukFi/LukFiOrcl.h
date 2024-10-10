/*--------------------------------------------------------------------------*/
/*---------------------------- File LukFiOrcl.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Definition of the LukFiOrcl class, which implements the FiOracle interface
 * and provides a set of 29 test problems for nonsmmoth unconstrained
 * optimization. The description of test problems is reported by
 * Ladislav Luksan and Jan Vlcek, it is available at
 *
 *     http://www.cs.cas.cz/~luksan/test.html
 *
 * and in
 *
 *     "Test problems for nonsmooth unconstrained optimization and
 *      linearity constrained optimization"
 *     Technical Report 798, Institute of Computer Science,
 *     Academy of Science of the Czech Republic
 *
 * \author  Annabella Astorino \n
 *          Istituto di Calcolo e Reti ad Alte Prestazioni - CNR \n
 *
 * \author  Antonio Frangioni \n
 *          Dipartimento di Informatica \n
 *          Universita' di Pisa \n
 *
 * \author  Manlio Gaudioso \n
 *          Dipartimento di Elettronica Informatica e Sistemistica \n
 *          Universita' della Calabria \n
 *
 * \author  Enrico Gorgone \n
 *          Dipartimento di Elettronica Informatica e Sistemistica \n
 *          Universita' della Calabria \n
 *
 * \copyright &copy; by Antonio Frangioni, Enrico Gorgone
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __LukFiOrcl
 #define __LukFiOrcl  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FiOracle.h"

/*--------------------------------------------------------------------------*/
/*------------------------ NAMESPACE and USINGS ----------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------------- CLASS LukFiOrcl ------------------------------*/
/*--------------------------------------------------------------------------*/

class LukFiOrcl : public FiOracle
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

 LukFiOrcl( std::istream *iStrm = nullptr );

/**< Constructor of the class. The parameter `iStrm', if provided, is taken
   as a pointer to a istream from which the algorithmic parameters for the
   subgradient algorithm are sequentially read in the following order. Each
   parameter must be placed at the beginning of a separate line, max 255
   characters long, with all the rest of the line up to the first newline
   character (apart from a separating whitespace) being available for
   comments. Any line whose first character is '#' and any blank line is
   ignored. If 0 is passed, the file ends before reaching a given parameter,
   or some parameter is in the wrong format, each non-specified parameter is
   given a default value, shown in [] below.

   NameF   [1]: the function that the FiOracle returns, chosen in the list
           of available 25 ones below

   NrVr    [2]: a few functions of the list, namely 26, 27 and 28, have a
           parametric number of variable, that is fixed to NrVr

   NrCmp   [2]: function 28 (MaxQR) has a parametric number of components,
           i.e., it is the sum of NrCmp functions

   seed    [1]: function 28 (MaxQR) is randomly generated, with this seed

   Here comes the list of functions:

   1 Rosenbrock. [ NonConvex ]. Optimal point \f$\lambda^* = (1, 1)\f$.
      Optimal value \f$f^* = 0\f$.
      \f[f(\lambda) = 100(\lambda_2-\lambda_1^2)^2+(1-\lambda_1)^2 \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   2 Crescent. [ NonConvex ]. Optimal point \f$\lambda^* = (0, 0)\f$.
      Optimal value \f$f^* = 0\f$.
      \f[f(\lambda) = \max\left\{ \lambda_1^2+(\lambda_2-1)^2 + \lambda_2 - 1\,
      ,\,-\lambda_1^2-(\lambda_2-1)^2 + \lambda_2 + 1 \right\} \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   3 CB2 (Charalambous/Bandler). [ Convex ].
      Optimal point \f$\lambda^* = (1.139286, 0.899365)\f$.
      Optimal value \f$f^* = 1.9522245\f$.
      \f[f(\lambda) = \max\left\{\lambda_1^2+\lambda_2^4\,,\, (2-\lambda_1)^2
      + (2-\lambda_2)^2\,,\,2e^{-\lambda_1+\lambda_2} \right\} \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   4 CB3 (Charalambous/Bandler). [ Convex ].
      Optimal point \f$\lambda^* = (1, 1)\f$. Optimal value \f$f^* = 2\f$.
      \f[f(\lambda) = \max\left\{\lambda_1^4+\lambda_2^2\,,\, (2-\lambda_1)^2
      + (2-\lambda_2)^2\,,\, 2e^{-\lambda_1+\lambda_2} \right\} \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   5 DEM (Demyanov/Malozemov). [ Convex ].
      Optimal point \f$\lambda^* = (0,-3)\f$. Optimal value \f$f^* = -3\f$
      \f[f(\lambda) = \max\left\{5\lambda_1+\lambda_2\,,\, -5\lambda_1+
      \lambda_2\,,\,\lambda_1^2+\lambda_2^2+4\lambda_2 \right\} \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   6 QL. [ Convex ].
      Optimal point \f$\lambda^* = (1.2, 2.4)\f$. Optimal value \f$f^* =
      7.2\f$.
      \f[f(\lambda) = \max\left\{\lambda_1^2+\lambda_2^2\,,\, \lambda_1^2+
      \lambda_2^2+ 10(-4\lambda_1-\lambda_2+4)\,,\,\lambda_1^2+\lambda_2^2
      + 10(-\lambda_1-2\lambda_2+6) \right\} \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   7 LQ. [ Convex ].
      Optimal point \f$\lambda^* = (\frac{1}{\sqrt 2}, \frac{1}{\sqrt 2})\f$.
      Optimal value \f$f^* = -\sqrt 2\f$.
      \f[f(\lambda) = \max\left\{-\lambda_1-\lambda_2\,,\, -\lambda_1-
      \lambda_2+(\lambda_1^2+\lambda_2^2-1)\right\} \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   8 Mifflin1. [ Convex ].
      Optimal point \f$\lambda^* = (1, 0)\f$. Optimal value \f$f^* = -1\f$.
      \f[f(\lambda) = -\lambda_1+20\max\left\{\lambda_1^2+\lambda_2^2-1 \,,
      \, 0\right\} \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   9 Mifflin2. [ NonConvex ].
      Optimal point \f$\lambda^* = (1, 0)\f$. Optimal value \f$f^* = -1\f$.
      \f[f(\lambda) = -\lambda_1+2(\lambda_1^2+\lambda_2^2-1)
      + 1.75\left|\lambda_1^2+\lambda_2^2-1\right| \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   10 Wolfe. [ NonConvex ]. Optimal value \f$f^* = -8\f$
      \f[ f(\lambda) = \left \{
      \begin{array}{ll} 5\sqrt{9\lambda_1^2+16\lambda_2^2} &
      \lambda_1\geq|\lambda_2| \\
      9\lambda_1+16|\lambda_2|& 0<\lambda_1<|\lambda_2|\\
      9\lambda_1+16|\lambda_2|-\lambda_1^9&\lambda_1\leq 0
      \end{array} \right. \f]
      See "Test problems for nonsmooth unconstrained optimization and
      linearity constrained optimization". Technical Report 798,
      Institute of Computer Science, Academy of Science of the Czech
      Republic, (2000)

   11 Rosen. [ Convex ].
      Optimal point \f$\lambda^* = (0, 1, 2, -1)\f$. Optimal value
      \f$f^* = -44\f$
      \f[f(\lambda) = \max\left\{ f_1(\lambda) \,,\,
      f_1(\lambda)+10f_2(\lambda)\,,\, f_1(\lambda)+10f_3(\lambda)\,,\,
      f_1(\lambda)+10f_4(\lambda)\right\}, \f]
      where
      \f[\left |
      \begin{array}{ll}
      f_1(\lambda)=& \lambda_1^2+\lambda_2^2+2\lambda_3^2+\lambda_4^2
      -5\lambda_1-5\lambda_2-21\lambda_3+7\lambda_4 \\
      f_2(\lambda)=& \lambda_1^2+\lambda_2^2+\lambda_3^2+\lambda_4^2
      +\lambda_1-\lambda_2+\lambda_3-\lambda_4-8\\
      f_3(\lambda)=& \lambda_1^2+2\lambda_2^2+\lambda_3^2+2\lambda_4^2
      -\lambda_1-\lambda_4-10\\
      f_4(\lambda)=& \lambda_1^2+\lambda_2^2+\lambda_3^2
      +2\lambda_1-\lambda_2-\lambda_4-5
      \end{array} \right. \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   12 Shor. [ Convex ].
      Optimal point \f$\lambda^* = (1.12434, 0.97945, 1.47770, 0.92023,
      1.12429)\f$. Optimal value \f$f^* = 22.60016\f$.
      \f[f(\lambda) = \max_{1\leq i \leq 10} \left\{ b_i \sum_{j=1}^5
      ( \lambda_{j} - a_{ij} )^2 \right\} \f]
      See "Methods of Descent for Nondifferentiable Optimization",
      K. Kiwiel, Lectures Notes in Mathematics, A. Dold and B. Eckmann Eds.
      Springer-Verlag, (1985)

   13  Maxquad (by Lemarechal). [ Convex ].
      Optimal point \f$\lambda^* = (-0.1263, -0.0346, -0.0067, -0.2668,
      0.0673, 0.2786, 0.0744, 0.1387, 0.0839, 0.0385)\f$.
      Optimal value \f$f^* = -0.8414084\f$.
      \f[f(\lambda) =  \max_{1\leq k \leq 5} \left( \lambda^TA_k \lambda -
      b_k ^T\lambda \right) \f]
      See "Nonsmooth Optimization", Proceedings of a IISA Workshop,
      C. Lemarechal and R. Mifflin, Eds. Vol. 3 (1977)

   14 Maxq. [ Convex ]. Optimal point \f$\lambda^*= 0\f$.
      Optimal value f* = 0.
      \f[f(\lambda) = \max_{1\leq i \leq 20} \lambda_i^2 \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   15 Maxl. [ Convex ]. Optimal point \f$\lambda^*= 0\f$.
      Optimal value \f$f^* = 0\f$.
      \f[f(\lambda) = \max_{1\leq i \leq 20} |\lambda_i| \f]
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   16 TR48 (by Goffin). [ Convex ].
      Optimal point \f$\lambda^*\f$ = (144, 257, 0, 483, 89, -165, -72, -252,
      -88, -178, 311, 126, 7, -135, 158, 209, 101, -92, 229, 80, 95, 71, -244,
      102, -12, 132, 337, 61, 104, 41, 261, 118, 99, -246, 156, -270, 330,
      -130, 952, -62, 161, 484, 122, 474, 1086, 861, -170, 206).
      Optimal value \f$f^* = -638565\f$.
      \f[f(\lambda) = \sum_{j=1}^{48} d_j\max_{1\leq i\leq 48}(\lambda_i-a_{ij})
      - \sum_{i=1}^{48} s_i\lambda_i \f]
      See "Nonsmooth Optimization", Proceedings of a IISA Workshop,
      C. Lemarechal and R. Mifflin, Eds. Vol. 3 (1977)

   17 Colville1. [ NonConvex ].
      Optimal value \f$f^* = -32.348679\f$.

   18 HS78. [ NonConvex ].
      Optimal value \f$f^* = -2.9197004\f$.

   19 El-Attar. [ NonConvex ].
      Optimal value \f$f^* = -0.5598131\f$. \n
      See "An algorithm for l1-norm minimization with application to nonlinear
      l1-approximation", R.A. El-Attar, M. Vidyasagar and S.R.K. Dutta,
      Siam Journal Numeircal Analysis, 16(1), pp. 70-86, (1979)

   20 Gill. [ NonConvex ].
      Optimal value \f$f^* = 9.7857721\f$.

   21 Steiner2. Optimal value \f$f^* = 16.703838\f$. \n

   22 Goffin. [ Convex ].
      Optimal value \f$f^* = 0\f$. \n
      See "Nonsmooth Optimization", M. Makela and P. Neittaanmaki, (1992)

   23 MXHILB. [ NonConvex ].
      Optimal value \f$f^* = 0\f$.

   24 L1HILB. [ NonConvex ].
      Optimal value \f$f^* = 0\f$. \n

   25 Shell Dual. [ NonConvex ].
      Optimal value \f$f^* = 32.348679\f$.

  26  smooth

  27  AbsVal

  28  MaxQR

  29  Lewis
  */

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the problem
    @{ */

   void GetInitialPoint( LMRow lmb );

/**< GetInitialPoint() writed into lmb[], which must be long at least
   GetNumVar(), the "standard starting point" of the function:

   @note Rosenbrock (NwNF==1). \f$\lambda = (-1.2, 1)\f$, for which
      \f$f = 24.2\f$.

   @note Crescent (NwNF==2). \f$\lambda = (-1.5, 2)\f$, for which
      \f$f = 4.25\f$.

   @note CB2 (NwNF==3). \f$\lambda = (1, -0.1)\f$, for which \f$f = 5.41\f$.

   @note CB3 (NwNF==4). \f$\lambda = (2, 2)\f$, for which \f$f = 20\f$.

   @note DEM (NwNF==5). \f$\lambda = (1, 1)\f$, for which \f$f = 6\f$.

   @note QL (NwNF==6). \f$\lambda = (-1, 5)\f$, for which \f$f = 56\f$.

   @note LQ (NwNF==7). \f$\lambda = (-0.5, -0.5)\f$, for which \f$f =1\f$.

   @note Mifflin1 (NwNF==8).
      \f$\lambda = (0.8, 0.6)\f$, for which \f$f =-0.8\f$.

   @note Mifflin2 (NwNF==9).
      \f$\lambda = (-1, -1)\f$, for which \f$f =4.75\f$.

   @note Wolfe. (NwNF==10). \f$\lambda = (3, 2)\f$, for which
      \f$f =60,207972894\f$.

   @note Rosen. (NwNF==11). \f$\lambda = (0, 0, 0, 0, 0)\f$, for which
      \f$f = 0\f$.

   @note Shor (NwNF==12).  \f$\lambda = (0, 0, 0, 0, 1)\f$, for which
      \f$f = 80\f$.

   @note Maxquad (NwNF==16). \f$\lambda = (0, 0, 0, 0, 0)\f$, for which
      \f$f = 0\f$
      [ \f$\lambda = (1, 1, 1, 1, 1)\f$, for which \f$f = 5337\f$ ].

   @note Maxq (NwNF==19). \f$\lambda_i =
      \left\{\begin{array}{ll} i& i=1,\ldots,10\\
      -i& i=11,\ldots,20\end{array} \right.\f$
      for which \f$f = 400\f$.

   @note Maxl (NwNF==20). \f$\lambda_i =
      \left\{\begin{array}{ll} i& i=1,\ldots,10\\
      -i& i=11,\ldots,20\end{array} \right.\f$
      for which \f$f = 20\f$.

   @note TR48 (NwNF==21). \f$\lambda\f$ = 0, for which \f$f = -464816\f$.
      [ \f$\lambda\f$ = (11.19, 127.2, -129.7, 344.5, -40.72, -295.3,
      -202.3, -382.3, -217.7, -307.7, 178.1, -4.36, -123.3, -265.3, 28.28,
      70.57, -31.81, -222.3, 96.19, -52.79, -34.71, -59.16, -373.7, -28.35,
      -141.7, 2.28, 198.5, -69.16, -26.35, -88.72, 130.8, -12.35, -30.7,
      -376.3, 23.18, -400.3, 197.1, -260.3, 813.5, -191.7, 31.29, 345.5,
      -7.72, 335.5, 947.5, 722.5, -300.3, 73.2) , for which
      \f$f = -638524.94\f$ ].

   ...

   17 Smooth. [ Convex ]. Somma dei mezzi quadrati:
      x0= 1 1 1 1 1 1  .

   18 AbsVal. [ Convex ]. Somma dei moduli:
      x0= 1 1 1 1 1 1 .. . ..

   19 Lewis
      x0 = -10 [1, 1,]
  */

/*--------------------------------------------------------------------------*/
/** Returns the name of the function (between 1 and 250. */

   Index GetName( void ) { return( NameF ); }

/**@ -----------------------------------------------------------------------*/
/*----------------------- METHODS FOR COMPUTING Fi() -----------------------*/
/*--------------------------------------------------------------------------*/

   virtual HpNum Fi( cIndex wFi = Inf< Index >() );

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR READING SUBGRADIENTS / CONSTRAINTS -------------*/
/*--------------------------------------------------------------------------*/

   virtual Index GetGi( SgRow SubG , cIndex_Set &SGBse ,
			cIndex Name = Inf< Index >() , cIndex strt = 0 ,
			Index stp = Inf< Index >() );

/*--------------------------------------------------------------------------*/

   virtual HpNum GetVal( cIndex Name = Inf< Index >() );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING OTHER RESULTS -------------------*/
/*--------------------------------------------------------------------------*/

   virtual HpNum GetLowerBound( cIndex wFi = Inf< Index >() );

/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

   ~LukFiOrcl();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED METHODS ---------------------------*/
/*--------------------------------------------------------------------------*/

   void SetDimension( void );

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
/*--------------------------------------------------------------------------*/

   Index NameF;       ///< name (number) of the function at point Lambda
   Index NrCmp;       ///< number of component functions
   Index seed;        ///< seed for random number generation

   Row bQR;
   Row aQR;
   LMMat cQR;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class LukFiOrcl )

/*--------------------------------------------------------------------------*/

};  // end( namespace NDO_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* LukFiOrcl.h */

/*--------------------------------------------------------------------------*/
/*------------------------- End File LukFiOrcl.h ---------------------------*/
/*--------------------------------------------------------------------------*/
