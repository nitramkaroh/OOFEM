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

#ifndef ogdenmaterial_h
#define ogdenmaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"

///@name Input fields for OgdenMaterial
//@{
#define _IFT_OgdenMaterial_Name "ogdenmat"
#define _IFT_OgdenMaterial_alpha "alpha"
#define _IFT_OgdenMaterial_mu "mu"
#define _IFT_OgdenMaterial_k "k"
//@}

namespace oofem {
/**
 * This class implements Compressible Ogden material.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 * 
 * References: R.W. Ogden: Non-Linear Elastic Deformations,
 * de Souza Neto, Peric, Owen: Computational Methods for Plasticity: Theory and Applications
 *
 * Free energy is considered as:
 * @f[
 * \rho_0 \psi = \sum_I^N \frac{mu_I}{alpha_I}(\lambda_1^{\alpha_I}+\lambda_2^{\alpha_I}+\lambda_3^{\alpha_I}-3) + \frac{1}{2} K[ln(J)]^2
 * @f]
 * @f$ \alpha_I @f$, @f$ mu_I @f$, and @f$K@f$ are material parameters.
 *
 * @f$ \lambda_1@f$, @f$ \lambda_2 @f$, and @f$ \lambda_3 @f$ are principal stretches
 *
 * Compressible Neo-Hookean model is obtained by setting @f$N = 1@f$, @f$\alpha_1 = 2@f$
*
 * Compressible Mooney-Rivlin model is obtained by setting @f$N = 2@f$, @f$\alpha_1 = 2@f$, and @f$\alpha_2 = -2@f$.
 */
class OgdenMaterial : public StructuralMaterial
{
protected:
    // Material parameters
    FloatArray alpha;
    FloatArray mu;
    double K;
    int N;


public:
    OgdenMaterial(int n, Domain *d);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }


    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep);


    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);

    
    void giveSecondPKStressVector_3d(FloatMatrix &answer, const FloatMatrix &C, const FloatArray &eVlas, const FloatMatrix &eVecs,TimeStep *tStep);
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual const char *giveInputRecordName() const { return _IFT_OgdenMaterial_Name; }
    virtual const char *giveClassName() const { return "OgdenMaterial"; }
};

} // end namespace oofem
#endif
