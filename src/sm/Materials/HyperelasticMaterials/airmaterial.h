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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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


#ifndef airmaterial_h
#define airmaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"



///@name Input fields for AirMaterial
//@{
#define _IFT_AirMaterial_Name "airmaterial"
#define _IFT_AirMaterial_p "p"
#define _IFT_AirMaterial_ltf "ltf"
#define _IFT_AirMaterial_eps "eps"
#define _IFT_AirMaterial_rmt "rmt"
//@}

namespace oofem {
/**
 * Free energy is considered as:
 * \f$[
 * \rho_0 \psi = pJ + eps d2W/dF2\f$
 *  This form of energy corresponds to an air material with regularization due to singularity of the shear stiffness.
 * @author Martin Horak, nitramkaroh@seznam.cz
 * @note Reference: article{@todo,
 * title={@todo},
 * author={@todo},
 * journal={@todo},
 * volume={46},
 * number={2},
 * pages={201--215},
 * year={1984},
 * publisher={Elsevier}
 * }
 *
 *
 */
class AirMaterial : public StructuralMaterial
{
protected:
    double p = 0.; ///< pressure
    int ltf  = 0.; ///< load time function
    double eps  = 0.; ///< regularization parameter
    StructuralMaterial *regularizationMaterial;

public:
    AirMaterial(int n, Domain *d);
    ~AirMaterial();

    IRResultType initializeFrom(InputRecord *ir);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
					    MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    
    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    
    void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);
    
    void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    
    MaterialStatus * CreateStatus(GaussPoint *gp) const;

    virtual const char *giveInputRecordName() const { return _IFT_AirMaterial_Name; }
    virtual const char *giveClassName() const { return "AirMaterial"; }
};
} // end namespace oofem
#endif
