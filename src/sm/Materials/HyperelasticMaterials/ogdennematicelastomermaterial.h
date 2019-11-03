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

#ifndef ogdennematicelastomermaterial_h
#define ogdennematicelastomermaterial_h

#include "../sm/Materials/HyperelasticMaterials/ogdenmaterial.h"

///@name Input fields for OgdenMaterial
//@{
#define _IFT_OgdenNematicElastomerMaterial_Name "ogdennematicmat"
#define _IFT_OgdenNematicElastomerMaterial_a "a"
#define _IFT_OgdenNematicElastomerMaterial_qce "qce"
//@}

namespace oofem {
/**
 * This class implements Compressible Ogden-type material for Nematic elastomer.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 * 
 * References: 
 * @article{agostiniani2012ogden,
 * title={Ogden-type energies for nematic elastomers},
 * author={Agostiniani, Virginia and DeSimone, Antonio},
 * journal={International Journal of Non-Linear Mechanics},
 * volume={47},
 * number={2},
 * pages={402--412},
 * year={2012},
 * publisher={Elsevier}
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
class OgdenNematicElastomerMaterial : public OgdenMaterial
{
protected:
    // Material parameters
    double a;
    int qcEnvelop;

public:
      OgdenNematicElastomerMaterial(int n, Domain *d);
      virtual IRResultType initializeFrom(InputRecord *ir);
      virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) override;
      virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep) override;

      virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

      virtual const char *giveInputRecordName() const { return _IFT_OgdenNematicElastomerMaterial_Name; }
      virtual const char *giveClassName() const { return "OgdenNematicElastomerMaterial"; }

 protected:
      void giveSecondPKStressVector_3d(FloatMatrix &answer, const FloatMatrix &C, const FloatArray &eVlas, const FloatMatrix &eVecs,TimeStep *tStep);
      void giveSecondPKStressVectorSmectic_3d(FloatMatrix &answer, const FloatMatrix &C, const FloatArray &eVlas, const FloatMatrix &eVecs,TimeStep *tStep);
      void give3dMaterialStiffnessMatrixSmectic_dSdE(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
      void give3dMaterialStiffnessMatrix_dSdE(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
      int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};
 


 

} // end namespace oofem
#endif
