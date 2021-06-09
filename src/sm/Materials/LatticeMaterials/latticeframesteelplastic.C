/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2019   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "LatticeFrameSteelPlastic.h"
#include "latticematstatus.h"
#include "latticestructuralmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "CrossSections/structuralcrosssection.h"
#include "engngm.h"
#include "mathfem.h"
#include "Elements/LatticeElements/latticestructuralelement.h"
#include "datastream.h"
#include "staggeredproblem.h"
#include "contextioerr.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LatticeFrameSteelPlastic);

// constructor which creates a dummy material without a status and without random extension interface
// LatticeFrameElastic :: LatticeFrameElastic(int n, Domain *d) :
//     LatticeStructuralMaterial(n, d)
// {}


bool
LatticeFrameSteelPlastic::hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void
LatticeFrameSteelPlastic::initializeFrom(InputRecord &ir)
{
    LatticeStructuralMaterial::initializeFrom(ir);

    //Young's modulus of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->e, _IFT_LatticeFrameSteelPlastic_e); // Macro

    //Poisson's ratio of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->nu, _IFT_LatticeFrameSteelPlastic_n); // Macro

    //N0
    IR_GIVE_FIELD(ir, this->N0, _IFT_LatticeFrameSteelPlastic_N0); // Macro

    //Mx0
    IR_GIVE_FIELD(ir, this->Mx0, _IFT_LatticeFrameSteelPlastic_Mx0); // Macro

    //My0
    IR_GIVE_FIELD(ir, this->My0, _IFT_LatticeFrameSteelPlastic_My0); // Macro

    //Mz0
    IR_GIVE_FIELD(ir, this->Mz0, _IFT_LatticeFrameSteelPlastic_Mz0); // Macro

    Tol== 1.e-6;
    IR_GIVE_FIELD(ir, this->Tol, _IFT_LatticeFrameSteelPlastic_tol); // Macro

    newtonIter = 100;
    IR_GIVE_FIELD(ir, this->newtonIter, _IFT_LatticeFrameSteelPlastic_iter); // Macro

    numberOfSubIncrements = 10;
    IR_GIVE_FIELD(ir, this->numberOfSubIncrements, _IFT_LatticeFrameSteelPlastic_sub); // Macro







    //Peter: You need to think which other material parameters you will need for your plastic model. I can think of the parameters that are used to normalise the forces in the yield condition (N0, Mx0, My0, Mz0) as well as a yieldtol. Maybe there are more. Add then the corresponding entries at the top of the h-file.
    
}

MaterialStatus *
LatticeFrameSteelPlastic::CreateStatus(GaussPoint *gp) const
{
    return new LatticeMaterialStatus(gp);
}

MaterialStatus *
LatticeFrameSteelPlastic::giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( !status ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status ) {
            gp->setMaterialStatus(status);
        }
    }

    return status;
}


FloatArrayF< 6 >
LatticeFrameSteelPlastic::giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const
//
// returns a FloatArray(6) of initial strain vector
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    double alpha = this->give(tAlpha, gp);


    return {
        alpha, 0., 0., 0., 0., 0.
    };
}


FloatArrayF< 6 >
LatticeFrameSteelPlastic::giveFrameForces3d(const FloatArrayF< 6 > &originalstrain,
                                       GaussPoint *gp,
                                       TimeStep *tStep)
{
    auto status = static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);
    /*Peter: Here you now need to work on the plastic return. You can check latticedamageplastic how it its done. Check in there perform plasticity return.
     */

    // First remove thermal strain from the original strain. 
    auto reducedStrain = originalStrain;
    auto thermalStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( thermalStrain.giveSize() ) {
        reducedStrain -= FloatArrayF< 6 >(thermalStrain);

    }

{
   

    /* Get plastic strain vector from status*/
    auto tempPlasticStrain = status->givePlasticLatticeStrain() [ { 0, 1, 2, 3 } ];

    FloatArrayF< 4 >tangent = { this->eNormalMean, this->alphaOne * this->eNormalMean, this->alphaOne * this->eNormalMean };
    /* Compute trial stress*/
    auto stress = mult(tangent, strain - tempPlasticStrain);

    //Introduce variables for subincrementation
    //Only _3dLattice is possible

    auto oldStrain = this->giveReducedStrain(gp, tStep) [ { 0, 1, 2, 3 } ];

    /* Compute yield value*/
    double yieldValue = computeYieldValue(stress, tempKappa, gp, tStep);
    int subIncrementCounter = 0;

    /* Check yield condition, i.e. if the yield value is less than the yield tolerance.
     * If yield condition is valid. Do perform regular return (closest point return)*/

    if ( yieldValue / pow(fcLocal, 2.) > yieldTol ) {
        // introduce a subincrementation flag
        int subIncrementFlag = 0;
        auto convergedStrain = oldStrain;
        auto tempStrain = strain;
        auto deltaStrain = strain - oldStrain;
        //To get into the loop
        returnResult = RR_NotConverged;
        while ( returnResult == RR_NotConverged || subIncrementFlag == 1 ) {
            stress = mult(tangent, tempStrain - tempPlasticStrain);
    /*Peter: Next you need to write the plasticity return. I would suggest to write it as it was done in latticeplasticitydamage in a separate function called performPlasticity return. The steps are 1) compute yield function based on subset of stresses. If larger than yield tolerance, perform return. You need to implement a lot of function such as computeYieldFunction, computeFVector (this stands for dFdSigma), computeMVector (this is DGDSigma, which is in your case not needed because you want to use associated flow, computeDMMatrix (which is for second derivative of F). The name of the functions is up to you. Please follow oofem programming rules, which means that you should give functions names which make sense as much as possible and use captital letters for every new word.  
     */

    //Peter: Here is an example of a call to the function performPlasticityReturn
    //    auto stress = this->performPlasticityReturn(gp, reducedStrain, tStep);

double
LatticeFrameSteelPlastic::performRegularReturn(FloatArrayF< 4 > &stress,
                                              double yieldValue,
                                              GaussPoint *gp,
                                              TimeStep *tStep) const
{
   // double fcLocal =  giveCompressiveStrength(gp, tStep);

    auto status = static_cast< LatticeFrameSteelPlasticStatus * >( this->giveStatus(gp) );

    double deltaLambda = 0.;

    auto trialStress = stress;
    auto tempStress = trialStress;

    //double trialShearStressNorm = norm(trialStress [ { 1, 2 } ]);

    //double tempShearStressNorm = trialShearStressNorm;

    //double thetaTrial = atan2(stress.at(3), stress.at(2) );

   
    //initialise unknowns
    FloatArrayF< 4 >unknowns;
    unknowns.at(1) = trialStress.at(1);
    //unknowns.at(2) = trialShearStressNorm;
    //unknowns.at(3) = tempKappa;
    unknowns.at(4) = 0.;

    // Look at the magnitudes of the residuals. You have to scale the yieldValue down.
    yieldValue = computeYieldValue(tempStress, tempKappa, gp, tStep);

    //initiate residuals
    FloatArrayF< 4 >residuals;
    residuals.at(4) = yieldValue;

    double normOfResiduals  = 1.; //just to get into the loop

    int iterationCount = 0;
    while ( normOfResiduals > yieldTol ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            returnResult = RR_NotConverged;
            return 0.;
        }

        //Normalize residuals. Think about it more.
        FloatArrayF< 4 >residualsNorm;
        residualsNorm.at(1) = residuals.at(1) / fcLocal;
        residualsNorm.at(2) = residuals.at(2) / fcLocal;
        residualsNorm.at(3) = residuals.at(3);
        residualsNorm.at(4) = residuals.at(4) / pow(fcLocal, 2.);

        normOfResiduals = norm(residualsNorm);

        //First check if return has failed
        if ( std::isnan(normOfResiduals) ) {
            returnResult = RR_NotConverged;
            return 0.;
        }

        if ( normOfResiduals > yieldTol ) {
            // Test to run newton iteration using inverse of Jacobian
            auto jacobian = computeJacobian(tempStress, tempKappa, deltaLambda, gp, tStep);

            auto solution = solve_check(jacobian, residuals);
            if ( solution.first ) {
                unknowns -= solution.second;
            } else {
                returnResult = RR_NotConverged;
                return kappa;
            }

            unknowns.at(2) = max(unknowns.at(2), 0.); //Keep rho greater than zero!
            unknowns.at(3) = max(unknowns.at(3), kappa); //Keep deltaKappa greater than zero!
            unknowns.at(4) = max(unknowns.at(4), 0.); //Keep deltaLambda greater than zero!

            /* Update increments final values and DeltaLambda*/
            tempStress.at(1) = unknowns.at(1);
            //tempShearStressNorm = unknowns.at(2);

            tempStress.at(2) = tempShearStressNorm * cos(thetaTrial);
            tempStress.at(3) = tempShearStressNorm * sin(thetaTrial);

            //tempKappa = unknowns.at(3);
            deltaLambda = unknowns.at(4);

            /* Compute the mVector holding the derivatives of the g function and the hardening function*/
            auto mVector = computeMVector(tempStress, tempKappa, gp, tStep);

            residuals.at(1) = tempStress.at(1) - trialStress.at(1) + this->eNormalMean * deltaLambda * mVector.at(1);
            residuals.at(2) = tempShearStressNorm - trialShearStressNorm + this->alphaOne * this->eNormalMean * deltaLambda * mVector.at(2);
            //residuals.at(3) = -tempKappa + kappa + deltaLambda * mVector.at(3);
            residuals.at(4) = computeYieldValue(tempStress, tempKappa, gp, tStep);
        }
    }

    returnResult = RR_Converged;

    stress = tempStress;

    return tempKappa;
}


    // auto stiffnessMatrix = LatticeFrameSteelPlastic::give3dFrameStiffnessMatrix(ElasticStiffness, gp, tStep);
    // auto stress = dot(stiffnessMatrix, strain);



    // printf("strain:\n");
    // strain.printYourself();
    // printf("stress:\n");
    // stress.printYourself();
    status->letTempLatticeStrainBe(strain);
    status->letTempLatticeStressBe(stress);

    return stress;
}


Interface *
LatticeFrameSteelPlastic::giveInterface(InterfaceType type)
{
    return nullptr;
}


FloatMatrixF< 6, 6 >
LatticeFrameSteelPlastic::give3dFrameStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime) const
{
    /*Peter: Gumaa, you need to enter here the part of the stiffness matrix which depend on the material parameters only. Later, you can then add the sectional parameters on the cross-section level*/
    static_cast< LatticeMaterialStatus * >( this->giveStatus(gp) );

    //Peter: Write here what shear modulus G is. It should be a function of E and nu.
    double g = this->e;

    //Peter: All the structural properties are read from the element
    const double area = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveArea();
    const double iy = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIy();
    const double iz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIz();
    const double ik = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIk();
    const double shearareay = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveShearAreaY();
    const double shearareaz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveShearAreaZ();

    //Peter: You need to put here the correct values. Please check this. 
    FloatArrayF< 6 >d = {
        this->e * area,
        g *shearareay,
        g *shearareaz,
        this->e * iy,
        this->e * iz,
        g *ik
    };
 
    return diag(d);
}
}
