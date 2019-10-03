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
 *               Copyright (C) 1993 - 2014   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef doublewellmaterial_h
#define doublewellmaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"

///@name Input fields for DoubleWellMaterial
//@{
#define _IFT_DoubleWellMaterial_Name "doublewellmaterial"
#define _IFT_DoubleWellMaterial_alpha "alpha"
#define _IFT_DoubleWellMaterial_eps "eps"
//@}

namespace oofem {
/**
 * This class implements Double Well Material.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 * 
 */
class DoubleWellMaterial : public StructuralMaterial
{
protected:
    // Material parameters
    double alpha;
    double eps;
    FloatMatrix tC1_0, tC2_0;

public:
    DoubleWellMaterial(int n, Domain *d);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }


    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep);


    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);

    virtual const char *giveInputRecordName() const { return _IFT_DoubleWellMaterial_Name; }
    virtual const char *giveClassName() const { return "DoubleWellMaterial"; }
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    
 protected:
    void compute_dC_dF(FloatMatrix &dCdF,const FloatArray &vF);
    FloatMatrix &givetC1(TimeStep *tStep);
    FloatMatrix &givetC2(TimeStep *tStep);

    
};

} // end namespace oofem
#endif
