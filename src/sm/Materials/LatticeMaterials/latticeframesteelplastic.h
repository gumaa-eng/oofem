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

#ifndef latticeframesteelplastic_h
#define latticeframesteelplastic_h

#include "latticestructuralmaterial.h"
#include "cltypes.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"
#include "latticematstatus.h"

///@name Input fields for LatticeFrameSteelPlastic
//@{
#define _IFT_LatticeFrameSteelPlastic_Name "LatticeFrameSteelPlastic"
#define _IFT_LatticeFrameSteelPlastic_talpha "talpha"
#define _IFT_LatticeFrameSteelPlastic_e "e"
#define _IFT_LatticeFrameSteelPlastic_n "n"
#define _IFT_LatticeFrameSteelPlastic_nx0 "nx0"
#define _IFT_LatticeFrameSteelPlastic_mx0 "mx0"
#define _IFT_LatticeFrameSteelPlastic_my0 "my0"
#define _IFT_LatticeFrameSteelPlastic_mz0 "mz0"
#define _IFT_LatticeFrameSteelPlastic_tol "tol"
#define _IFT_LatticeFrameSteelPlastic_iter "iter"
#define _IFT_LatticeFrameSteelPlastic_sub "sub"
//@}

namespace oofem {
/**
 * This class implements a local random linear elastic model for lattice elements.
 */
class LatticeFrameSteelPlastic : public LatticeStructuralMaterial

{

protected:
    ///Normal modulus
    double e;

    ///Ratio of shear and normal modulus
    double nu;

   ///nx0
    double nx0;

   ///mx0
    double mx0;

   ///my0
    double my0;

   ///mz0
    double mz0;

   ///tol
    double yieldTol;

   ///iter
    double newtonIter;

   ///sub
    double numberOfSubIncrements;

    enum LatticePlasticityDamage_ReturnResult { RR_NotConverged, RR_Converged };
    mutable LatticePlasticityDamage_ReturnResult returnResult = RR_NotConverged; /// FIXME: This must be removed. Not thread safe. Shouldn't be stored at all.

    double initialYieldStress = 0.;

    //

public:
    LatticeFrameSteelPlastic(int n, Domain *d) : LatticeStructuralMaterial(n, d) { };

    FloatArrayF< 4 >computeFVector(const FloatArrayF< 4 > &sigma, GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 4, 4 >computeDMMatrix(const FloatArrayF< 4 > &sigma, GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const override;

    FloatArrayF< 6 >giveReducedLatticeStrain(GaussPoint *gp, TimeStep *tStep) const;

    FloatArrayF< 6 >performPlasticityReturn(GaussPoint *gp, const FloatArrayF< 6 > &reducedStrain, TimeStep *tStep) const;

    double performRegularReturn(FloatArrayF< 4 > &stress, double yieldValue, GaussPoint *gp, TimeStep *tStep) const;

    double computeYieldValue(const FloatArrayF< 4 > &sigma, GaussPoint *gp, TimeStep *tStep) const;

    const char *giveInputRecordName() const override { return _IFT_LatticeFrameSteelPlastic_Name; }

    const char *giveClassName() const override { return "LatticeFrameSteelPlastic"; }

    void initializeFrom(InputRecord &ir) override;


    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }


    FloatMatrixF< 6, 6 >give3dFrameStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;


    bool hasMaterialModeCapability(MaterialMode mode) const override;


    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

protected:
};
} // end namespace oofem

#endif
