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

#include "latticeframeplasticsteel.h"
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
REGISTER_Material(LatticeFramePlasticsteel);

// constructor which creates a dummy material without a status and without random extension interface
// LatticeFrameElastic :: LatticeFrameElastic(int n, Domain *d) :
//     LatticeStructuralMaterial(n, d)
// {}


bool
LatticeFramePlasticsteel::hasMaterialModeCapability(MaterialMode mode) const
{
    return ( mode == _3dLattice );
}


void
LatticeFramePlasticsteel::initializeFrom(InputRecord &ir)
{
    LatticeStructuralMaterial::initializeFrom(ir);

    //Young's modulus of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->e, _IFT_LatticeFramePlasticsteel_e); // Macro

    //Poisson's ratio of the material that the beam element is made of
    IR_GIVE_FIELD(ir, this->nu, _IFT_LatticeFramePlasticsteel_n); // Macro

    //Peter: You need to think which other material parameters you will need for your plastic model. I can think of the parameters that are used to normalise the forces in the yield condition (N0, Mx0, My0, Mz0) as well as a yieldtol. Maybe there are more. Add then the corresponding entries at the top of the h-file.
    
}

MaterialStatus *
LatticeFramePlasticsteel::CreateStatus(GaussPoint *gp) const
{
    return new LatticeMaterialStatus(gp);
}

MaterialStatus *
LatticeFramePlasticsteel::giveStatus(GaussPoint *gp) const
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
LatticeFramePlasticsteel::giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const
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
LatticeFramePlasticsteel::giveFrameForces3d(const FloatArrayF< 6 > &originalstrain,
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

    /*Peter: Next you need to write the plasticity return. I would suggest to write it as it was done in latticeplasticitydamage in a separate function called performPlasticity return. The steps are 1) compute yield function based on subset of stresses. If larger than yield tolerance, perform return. You need to implement a lot of function such as computeYieldFunction, computeFVector (this stands for dFdSigma), computeMVector (this is DGDSigma, which is in your case not needed because you want to use associated flow, computeDMMatrix (which is for second derivative of F). The name of the functions is up to you. Please follow oofem programming rules, which means that you should give functions names which make sense as much as possible and use captital letters for every new word.  
     */

    //Peter: Here is an example of a call to the function performPlasticityReturn
    
    //    auto stress = this->performPlasticityReturn(gp, reducedStrain, tStep);
    

    // auto stiffnessMatrix = LatticeFramePlasticsteel::give3dFrameStiffnessMatrix(ElasticStiffness, gp, tStep);
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
LatticeFramePlasticsteel::giveInterface(InterfaceType type)
{
    return nullptr;
}


FloatMatrixF< 6, 6 >
LatticeFramePlasticsteel::give3dFrameStiffnessMatrix(MatResponseMode rmode, GaussPoint *gp, TimeStep *atTime) const
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
