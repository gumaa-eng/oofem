/* $Header: /home/cvs/bp/oofem/oofemlib/src/nrsolver.h,v 1.9.4.1 2004/04/05 15:19:43 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

  //   ***********************************
  //   *** CLASS NEWTON RAPHSON SOLVER ***
  //   ***********************************

 
#ifndef nrsolver_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "flotarry.h"
#include "cltypes.h"
#include "linesearch.h"
#ifdef __PETSC_MODULE
#ifndef __MAKEDEPEND
#include "petscksp.h"
#endif
#endif


  class Domain; class EngngModel;

/**
   This class implements the class NumericalMethod instance Newton-Raphson Method
   for solving non-linear problems. 
   Implements direct displacement controll without requiring BC applied. 
   EXPERIMENTAL!
   This is achieved by adding a large number alpha to the corresponding 
   diagonal member of K and replacing the right-hand side by alpha*prescribed_value.
   If alpha is very much larger than other stiffness corefficients then this alteration 
   effectively replaces the corresponding equation by 
   alpha*unknown_value = alpha*prescribed_value
   that is, the reqyuired condition, but the whole system remains symmetric and minimal
   changes are necessary in the computational sequence.
   The above artifice has been introduced by Payne and Irons.
*/
class NRSolver : public SparseNonLinearSystemNM
{
  /*
    DESCRIPTION :
    Perform solution of non-linear problem.

    =======>   This method uses Modified Newton Raphson iteration scheme  <======
   
    If we solve non-linear static we can interprete symbols as follows:

    Kt     - tangential stiffness
    deltaR - increment of displacements
    g      - vector of unballanced forces (at the end should be zero one)
    R      - Load vector (Incremental)
    R0     - Initial Load vector
    RT     - TotalLoadVector
    r      - total displacement vector
    F(r)   - Nodal representation of (real) internal forces.
    NR_Mode- variable controlling the mode of NRM (ModifiedNR, Full NRM (stifnees update after each iteration), 
    Modified Accelerated NRM (we perform iteration with stiffness matrix updated only after calm_MANRMSteps)
    calm_NR_OldMode - variable containing the old mode of NRM, which will be restored after 
    calm_NR_ModeTick iterations.
    calm_NR_ModeTick - see calm_NR_OldMode.
    calm_MANRMSteps - if calm_NR_Mode == calm_accelNRM, it specifies, that new updated
    stiffness matrix is assembled after calm_MANRMSteps.

    The load level and corresponding load vector is determined using intrinsic time, which is generated by
    nonlinear static model.


    TASKS :

    - solving problem 
    solveYourselfAt.
    - returning results (increment of displacement,
    reached level of loading and so on)

    Variable description  :

    K(N,N)    = STIFFNESS MATRIX (ASSUMED POZITIVE DEFINITE)        *
    deltaR(N) = ITERATIVE INCREMENT OF DISPLACEMENT                 *
    R         = LOAD VECTOR (Incr.)d                                *
    R0        = Initial Load Vector                                 *
    RT        = Total Load Vector                                   *
    DeltaR    = CURRENT TOTAL INCREMENT                             *
    F         = NODAL REPRESENTATION OF (REAL) INTERNAL FORCES      *

    RTOL      = CONVERGENCE TOLERANCE                               *

    OUTPUT : (after call solveYourselfAt)
    K(N,N)    = DIAGONALIZED STIFFNESS MATRIX                       *
    DeltaR    = REACHED DISPLACEMENT INCREMENT                      *
    nite      = NUMBER OF ITERATIONS REQUIRED TO FULLFIL BALANCE    *
    status    = NM_status with flags set to reached state (see cltypes.h) *
 
  */
 private:

  enum    nrsolver_ModeType {nrsolverModifiedNRM, nrsolverFullNRM, nrsolverAccelNRM};


  int            nite,nsmax;
  double         rtol, deltaL;
  double         minStepLength;
  int            solved;
  nrsolver_ModeType NR_Mode, NR_OldMode;
  int            NR_ModeTick;
  int            MANRMSteps;
  /// linear system solver
  SparseLinearSystemNM* linSolver;
  /// linear system solver ID
  LinSystSolverType solverType;
  /// sparse matrx version, used to controll constrains application to stiffness
  SparseMtrx::SparseMtrxVersionType smConstraintVersion;
  /// number of prescribed displacement
  int numberOfPrescribedDofs;
  /** flag indicating that some dofs are controlled under displa controoll.
      In parallel mode, numberOfPrescribedDofs is local (related to specific partition)
      so its nonzero value does not mean that there are no prescribed dofs on
      other partitions */
  bool prescribedDofsFlag;
  
  /// array of pairs identifying prescribed dofs (node, dof)
  IntArray prescribedDofs;
  /// array of prescribed values
  FloatArray prescribedDofsValues;
  /// load Time Function of prescribed values
  int prescribedDisplacementLTF;
  /// array of prescribed equations
  IntArray prescribedEqs;
  /// flag indicating that prescribedEqs were initialized
  bool prescribedEqsInitFlag;
  /// computed reactions. They are stored in order to print them in printState method.
  FloatArray lastReactions;
  /// flag indicating whether to use line-search
  int lsFlag;
  /// lineseach solver
  LineSearchNM* linesearchSolver;
#ifdef __PETSC_MODULE
  IS prescribedEgsIS;
  bool prescribedEgsIS_defined;
#endif


  public :
    NRSolver (int i, Domain* d,EngngModel* m, EquationID ut);
  // constructor
  ~NRSolver () ;              // destructor

  // solving
  /**
     Solves the given sparse linear system of equations g(x,l)=l-F(x); dx=K^{-1}g+ dl K^{-1}R.
     Total load vector not passed, it is defined as l*R+R0, where l is scale factor
     @param K coefficient matrix (K = dF/dx; stiffness matrix)
     @param R  incremental Rhs (incremental load)
     @param R0 initial Rhs (initial load)
     @param Rr linearization of K*rri, where rri is increment of prescribed displacements
     @param r  total solution (total displacement)
     @param dr increment of solution (incremental displacaments)
     @param l  Rhs scale factor (load level)
     @param rtol prescribed tolerance (g residual and iterative r change;)
     @param rlm - reference load mode
     @param F  InternalRhs (real internal forces)
     @return NM_Status value
  */
  virtual NM_Status solve (SparseMtrx* k, FloatArray* R, FloatArray* R0,
                           FloatArray* Rr, FloatArray* r, FloatArray* dr, FloatArray* F,
                           double& l, double rtol, referenceLoadInputModeType rlm,
                           int& nite, TimeStep*) ;

  virtual double giveCurrentStepLength() {return deltaL;}
  virtual void   setStepLength(double l) {deltaL = l;}
  /**
     Prints status mesage of receiver to output stream.
     Prints the message corresponding to last solve.
  */
  virtual void   printState (FILE* outputStream);


  // management  components
  IRResultType initializeFrom (InputRecord* ir);
  /** Stores receiver state to output stream. 
      Receiver should write class-id first in order to allow test
      whether correct data are then restored.
      @param stream output stream 
      @param mode determines ammount of info required in stream (state, definition,...)
      @param obj special parameter, used only to send particular integration
      point to material class version of this method. Except this 
      case, obj parameter is always NULL pointer.*/
  contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
  /** Restores the receiver state previously written in stream.
      @see saveContext member function.*/
  contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);

  // identification 
  const char*  giveClassName () const { return "NRSolver" ;}
  classType giveClassID () const { return NRSolverClass ;}
  /// sets associated Domain 
  virtual void         setDomain (Domain* d) {this->domain = d; if (linSolver) linSolver->setDomain(d);}
  /// This method clears receiver cached data dependent on topology, when it changes.
  virtual void reinitialize () {if (linSolver) linSolver->reinitialize();}
 protected:

  SparseLinearSystemNM* giveLinearSolver() ;
  LineSearchNM* giveLineSearchSolver() ;
  void initPrescribedEqs ();
  //void computeBCLoadVector (FloatArray& answer, SparseMtrx* k, TimeStep* atTime);
  void applyConstraintsToStiffness (SparseMtrx* k);
  void applyConstraintsToLoadIncrement (int nite, const SparseMtrx* k, FloatArray& R, 
                                        referenceLoadInputModeType rlm, TimeStep* atTime);

};

#define nrsolver_h
#endif









