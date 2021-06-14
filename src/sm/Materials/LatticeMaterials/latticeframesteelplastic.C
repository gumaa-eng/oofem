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

#include "latticeframesteelplastic.h"
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

  //gumaa
void
LatticeFrameSteelPlastic::initializeFrom(InputRecord &ir)
{
    LatticeStructuralMaterial::initializeFrom(ir);

    //Young's modulus of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->e, _IFT_LatticeFrameSteelPlastic_e); // Macro

    //Poisson's ratio of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->nu, _IFT_LatticeFrameSteelPlastic_n); // Macro

    //N0
    IR_GIVE_FIELD(ir, this->n0, _IFT_LatticeFrameSteelPlastic_n0); // Macro

    //Mx0
    IR_GIVE_FIELD(ir, this->mx0, _IFT_LatticeFrameSteelPlastic_mx0); // Macro

    //My0
    IR_GIVE_FIELD(ir, this->my0, _IFT_LatticeFrameSteelPlastic_my0); // Macro

    //Mz0
    IR_GIVE_FIELD(ir, this->mz0, _IFT_LatticeFrameSteelPlastic_mz0); // Macro

    tol= 1.e-6;
    IR_GIVE_FIELD(ir, this->tol, _IFT_LatticeFrameSteelPlastic_tol); // Macro

    iter = 100;
    IR_GIVE_FIELD(ir, this->iter, _IFT_LatticeFrameSteelPlastic_iter); // Macro

    sub = 10;
    IR_GIVE_FIELD(ir, this->sub, _IFT_LatticeFrameSteelPlastic_sub); // Macro







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

  //gumaa
FloatArrayF< 4 >
LatticeFrameSteelPlastic::computeFVector(const FloatArrayF< 4 > &stress, GaussPoint *gp, TimeStep *tStep) const

{
   
    double nx = nx;
    double mx = mx;
    double my = my;
    double mz = mz;

    FloatArrayF< 4 >f;
    {
        f.at(1) = 2.*(nx/this->n0)*(1/this->n0);
        f.at(2) = 2.*(mx/this->mx0)*(1/this->mx0);
        f.at(3) = 2.*(my/this->my0)*(1/this->my0);
        f.at(4) = 2.*(mz/this->mz0)*(1/this->mz0);
    
    }

    return f;
}

  //gumaa
FloatMatrixF< 4, 4 >
LatticeFrameSteelPlastic::computeDMMatrix(const FloatArrayF< 4 > &stress, GaussPoint *gp, TimeStep *tStep) const

{
    double nx = nx;
    double mx = mx;
    double my = my;
    double mz = mz;

    FloatMatrixF< 4, 4 >dm;
    { //main ellipse
        //Derivatives of dGDSig
        dm.at(1, 1) = 2/pow(this->n0, 2);
        dm.at(1, 2) = 0;
        dm.at(1, 3) = 0;
        dm.at(1, 4) = 0;

        //Derivatives of dGDTau
        dm.at(2, 1) = 0;
        dm.at(2, 2) = 2/pow(this->mx0, 2);
        dm.at(2, 3) = 0;
        dm.at(2, 4) = 0;

        //Derivates of evolution law
        dm.at(3, 1) = 0;
        dm.at(3, 2) = 0;
        dm.at(3, 3) = 2/pow(this->my0, 2);
        dm.at(3, 4) = 0;

        //Derivates of evolution law
        dm.at(4, 1) = 0;
        dm.at(4, 2) = 0;
        dm.at(4, 3) = 0;
        dm.at(4, 3) = 2/pow(this->mz0, 2);
    } 
    return dm;
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

  //gumaa
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
    auto reducedStrain = originalstrain;
    auto thermalStrain = this->computeStressIndependentStrainVector(gp, tStep, VM_Total);
    if ( thermalStrain.giveSize() ) {
        reducedStrain -= FloatArrayF< 6 >(thermalStrain);
       }
    //gumaa
    //Subset of reduced strain.
    //shear components are not used for plasticity return
    auto strain = reducedStrain [ { 0, 3, 4, 5 } ];
    const double area = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveArea();
    const double iy = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIy();
    const double iz = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIz();
    const double ik = ( static_cast< LatticeStructuralElement * >( gp->giveElement() ) )->giveIk();


    /* Get plastic strain vector from status*/
    auto tempPlasticStrain = status->givePlasticLatticeStrain() [ { 0, 3, 4, 5 } ];
    FloatArrayF< 4 >tangent = { this->e*area , this->e*ik, this->e*iy, this->e*iz };

    /* Compute trial stress*/
    auto stress = mult(tangent, strain - tempPlasticStrain);

    //Introduce variables for subincrementation
    //Only _3dLattice is possible

    auto oldStrain = this->giveReducedStrain(gp, tStep) [ { 0, 3, 4, 5 } ];

    /* Compute yield value*/
    double yieldValue = computeYieldValue(stress, gp, tStep);
    int subIncrementCounter = 0;

    /* Check yield condition, i.e. if the yield value is less than the yield tolerance.
     * If yield condition is valid. Do perform regular return (closest point return)*/

    if ( yieldValue > tol ) {
        // introduce a subincrementation flag
        int subIncrementFlag = 0;
        auto convergedStrain = oldStrain;
        auto tempStrain = strain;
        auto deltaStrain = strain - oldStrain;
        //To get into the loop
        returnResult = RR_NotConverged;
        while ( returnResult == RR_NotConverged || subIncrementFlag == 1 ) {
            stress = mult(tangent, tempStrain - tempPlasticStrain);

            if ( returnResult == RR_NotConverged ) {
                subIncrementCounter++;
                if ( sub > numberOfSubIncrements ) {
                    OOFEM_LOG_INFO("Unstable element %d \n", gp->giveElement()->giveGlobalNumber() );
                    OOFEM_LOG_INFO("Yield value %e \n", yieldValue);
                    OOFEM_LOG_INFO("ConvergedStrain value %e %e %e\n", convergedStrain.at(1), convergedStrain.at(2), convergedStrain.at(3), convergedStrain.at(4) );
                    OOFEM_LOG_INFO("tempStrain value %e %e %e %e\n", tempStrain.at(1), tempStrain.at(2), tempStrain.at(3), tempStrain.at(4) );
                    OOFEM_LOG_INFO("deltaStrain value %e %e %e %e\n", deltaStrain.at(1), deltaStrain.at(2), deltaStrain.at(3), deltaStrain.at(4) );
                    OOFEM_LOG_INFO("targetstrain value %e %e %e %e\n", strain.at(1), strain.at(2), strain.at(3), strain.at(4) );

                    OOFEM_ERROR("LatticeFrameSteelPlastic :: performPlasticityReturn - Could not reach convergence with small deltaStrain, giving up.");
                }
                printf("subincrementation required\n");
                subIncrementFlag = 1;
                deltaStrain *= 0.5;
                tempStrain = convergedStrain + deltaStrain;
            } else if ( returnResult == RR_Converged && subIncrementFlag == 1 ) {
                tempPlasticStrain.at(1) = tempStrain.at(1) - stress.at(1) / (area * this->e);
                tempPlasticStrain.at(2) = tempStrain.at(2) - stress.at(2) / ( ik * this->e );
                tempPlasticStrain.at(3) = tempStrain.at(3) - stress.at(3) / ( iy * this->e );
                tempPlasticStrain.at(4) = tempStrain.at(4) - stress.at(4) / ( iz * this->e );

                status->letTempPlasticLatticeStrainBe(assemble< 6 >(tempPlasticStrain, { 0, 3, 4, 5 }) );

                subIncrementFlag = 0;
                returnResult = RR_NotConverged;
                convergedStrain = tempStrain;
                deltaStrain = strain - convergedStrain;
                tempStrain = strain;
                subIncrementCounter = 0;
            } 
        }
    }


    tempPlasticStrain.at(1) = strain.at(1) - stress.at(1) / (area * this->e);
    tempPlasticStrain.at(2) = strain.at(2) - stress.at(2) / ( ik * this->e );
    tempPlasticStrain.at(3) = strain.at(3) - stress.at(3) / ( iy * this->e );
    tempPlasticStrain.at(4) = strain.at(4) - stress.at(4) / ( iz * this->e );

    status->letTempPlasticLatticeStrainBe(assemble< 6 >(tempPlasticStrain, { 0, 3, 4, 5  }) );

    return answer;
}
    /*Peter: Next you need to write the plasticity return. I would suggest to write it as it was done in latticeplasticitydamage in a separate function called performPlasticity return. The steps are 1) compute yield function based on subset of stresses. If larger than yield tolerance, perform return. You need to implement a lot of function such as computeYieldFunction, computeFVector (this stands for dFdSigma), computeMVector (this is DGDSigma, which is in your case not needed because you want to use associated flow, computeDMMatrix (which is for second derivative of F). The name of the functions is up to you. Please follow oofem programming rules, which means that you should give functions names which make sense as much as possible and use captital letters for every new word.  
     */

    //Peter: Here is an example of a call to the function performPlasticityReturn
    //    auto stress = this->performPlasticityReturn(gp, reducedStrain, tStep);
  //gumaa
double
LatticeFrameSteelPlastic::performRegularReturn(FloatArrayF< 4 > &stress,
                                              double yieldValue,
                                              GaussPoint *gp,
                                              TimeStep *tStep) const
{
    

    double deltaLambda = 0.;

    auto trialStress = stress;
    auto tempStress = trialStress;

    double trialStress = trialStress [ { 1, 2, 3, 4 } ]);

    double tempStress = trialStress;


    //initialise unknowns
    FloatArrayF< 4 >unknowns;
    unknowns.at(1) =           .at(1);
    unknowns.at(2) =           .at(2);
    unknowns.at(3) =           .at(3);
    unknowns.at(4) =           .at(4);
    unknowns.at(5) = 0;

    // Look at the magnitudes of the residuals. You have to scale the yieldValue down.
    yieldValue = computeYieldValue(tempStress, tempKappa, gp, tStep);

    //initiate residuals
    FloatArrayF< 5 >residuals;
    residuals.at(5) = yieldValue;

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
        residualsNorm.at(1) = residuals.at(1);
        residualsNorm.at(2) = residuals.at(2);
        residualsNorm.at(3) = residuals.at(3);
        residualsNorm.at(4) = residuals.at(4);
        residualsNorm.at(5) = residuals.at(5);


        normOfResiduals = norm(residualsNorm);

        //First check if return has failed
        if ( std::isnan(normOfResiduals) ) {
            returnResult = RR_NotConverged;
            return 0.;
        }

       
            /* Compute the mVector holding the derivatives of the g function and the hardening function*/
            auto mVector = computeMVector(tempStress, tempKappa, gp, tStep);

            residuals.at(1) = tempStress.at(1) - trialStress.at(1) + this->e * deltaLambda * fVector.at(1);
            residuals.at(2) = tempStress.at(2) - trialStress.at(2) + ik*this->e * deltaLambda * fVector.at(2);
            residuals.at(3) = tempStress.at(3) - trialStress.at(3) + iy*this->e * deltaLambda * fVector.at(3);
            residuals.at(4) = tempStress.at(4) - trialStress.at(4) + iz*this->e * deltaLambda * fVector.at(4);
            residuals.at(5) = computeYieldValue(tempStress, gp, tStep);

        }
    }

    returnResult = RR_Converged;

    stress = tempStress;

    
}


    // auto stiffnessMatrix = LatticeFrameSteelPlastic::give3dFrameStiffnessMatrix(ElasticStiffness, gp, tStep);
    // auto stress = dot(stiffnessMatrix, strain);



    // printf("strain:\n");
    // strain.printYourself();
    // printf("stress:\n");
    // stress.printYourself();

    status->letTempLatticeStrainBe(strain);
    status->letTempLatticeStressBe(stress);


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

