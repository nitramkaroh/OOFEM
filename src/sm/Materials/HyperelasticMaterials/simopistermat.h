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


#ifndef simopistermaterial_h
#define simopistermaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"



///@name Input fields for SimoPisterMaterial
//@{
#define _IFT_SimoPisterMaterial_Name "simopistermat"
#define _IFT_SimoPisterMaterial_mu "mu"
#define _IFT_SimoPisterMaterial_lambda "lambda"
//@}

namespace oofem {
/**
 * Free energy is considered as:
 * \f$[
 * \rho_0 \psi = U(J) + G( 0.5 * I_1 - ln J) ]\f$
 *  This form of energy corresponds to a neo-Hookean material which is extended to the compressible * range by adding an extra function depending on J.
 * @author Martin Horak, nitramkaroh@seznam.cz
 * @note Reference: article{simo1984remarks,
 * title={Remarks on rate constitutive equations for finite deformation problems: computational implications},
 * author={Simo, Juan C and Pister, Karl S},
 * journal={Computer Methods in Applied Mechanics and Engineering},
 * volume={46},
 * number={2},
 * pages={201--215},
 * year={1984},
 * publisher={Elsevier}
 * }
 *
 *
 */
class SimoPisterMaterial : public StructuralMaterial
{
protected:
    double lambda = 0.; ///< First Lame parameter
    double mu = 0.; ///< Second Lame parameter

public:
    SimoPisterMaterial(int n, Domain *d);

    IRResultType initializeFrom(InputRecord *ir);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
					    MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    
    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    
    void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);
    
    void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    
    MaterialStatus * CreateStatus(GaussPoint *gp) const;

    virtual const char *giveInputRecordName() const { return _IFT_SimoPisterMaterial_Name; }
    virtual const char *giveClassName() const { return "SimoPisterMaterial"; }
};
} // end namespace oofem
#endif
