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

#ifndef LatticeFrameSteelPlastic_h
#define LatticeFrameSteelPlastic_h

#include "latticestructuralmaterial.h"
#include "cltypes.h"
#include "randommaterialext.h"
#include "strainvector.h"
#include "stressvector.h"
#include "latticematstatus.h"

///@name Input fields for LatticeFrameSteelPlastic
//@{
#define _IFT_LatticeFrameSteelPlastic_Name "latticeframesteelplastic"
#define _IFT_LatticeFrameSteelPlastic_talpha "talpha"
#define _IFT_LatticeFrameSteelPlastic_e "e"
#define _IFT_LatticeFrameSteelPlastic_n "n"
#define _IFT_LatticeFrameSteelPlastic_n0 "n0"
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

   ///n0
    double n0;

   ///mx0
    double mx0;

   ///my0
    double my0;

   ///mz0
    double mz0;

   ///tol
    double tol;

   ///iter
    double iter;

   ///sub
    double sub;




public:
    LatticeFrameSteelPlastic(int n, Domain *d) : LatticeStructuralMaterial(n, d) { };

    FloatArrayF< 4 >computeFVector(const FloatArrayF< 4 > &sigma,
                                   GaussPoint *gp, TimeStep *tStep) const;

    FloatMatrixF< 4, 4 >computeDMMatrix(const FloatArrayF< 4 > &sigma, 
                                        GaussPoint *gp, TimeStep *tStep) const;



    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp,  TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_LatticeFrameSteelPlastic_Name; }

    const char *giveClassName() const override { return "LatticeFrameSteelPlastic"; }

    void initializeFrom(InputRecord &ir) override;


    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    FloatArrayF< 6 >giveFrameForces3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) override;

    FloatMatrixF< 6, 6 >give3dFrameStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const override;


    bool hasMaterialModeCapability(MaterialMode mode) const override;

    Interface *giveInterface(InterfaceType) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

protected:
};
} // end namespace oofem

#endif
