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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

#ifndef micropolarmaterial_elastic_h
#define micropolarmaterial_elastic_h


#include "../sm/Materials/Micromorphic/micromorphicmaterialextensioninterface.h"
#include "../sm/Materials/Micromorphic/micromorphicms.h"
#include "../sm/Materials/isolinearelasticmaterial.h"
#include "../sm/Materials/structuralmaterial.h"
#include "cltypes.h"


///@name Input fields for MicromorphLEmat
//@{
#define _IFT_MicropolarMaterial_Chiral_Name "micropolarmaterial_chiral"

#define _IFT_MicropolarMaterial_Chiral_lambda "lambda"
#define _IFT_MicropolarMaterial_Chiral_mu "mu"
#define _IFT_MicropolarMaterial_Chiral_mu_c "mu_c"

#define _IFT_MicropolarMaterial_Chiral_alpha "alpha"
#define _IFT_MicropolarMaterial_Chiral_beta "beta"
#define _IFT_MicropolarMaterial_Chiral_gamma "gamma"
#define _IFT_MicropolarMaterial_Chiral_C1 "c1"
#define _IFT_MicropolarMaterial_Chiral_C2 "c2"
#define _IFT_MicropolarMaterial_Chiral_C3 "c3"

//@}

namespace oofem {


/**
 * MicromorphicLinearElasticMaterial
 */
class MicropolarMaterial_Chiral : public IsotropicLinearElasticMaterial, MicromorphicMaterialExtensionInterface
{
protected:

  double lambda;
  double mu;
  double mu_c;

  double alpha;
  double beta;
  double gamma;


  double C1;
  double C2;
  double C3;

public:
    MicropolarMaterial_Chiral(int n, Domain * d);
    virtual ~MicropolarMaterial_Chiral();

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_MicropolarMaterial_Chiral_Name; }
    virtual const char *giveClassName() const { return "MicropolarMaterial_Chiral"; }

    virtual IRResultType initializeFrom(InputRecord *ir) override ;

    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == MicromorphicMaterialExtensionInterfaceType ) {
            return static_cast< MicromorphicMaterialExtensionInterface * >(this);
        } else {
            return NULL;
        }
    }

    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    virtual void giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override ;
    virtual void giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override ;
    virtual void giveMicromorphicMatrix_dSigdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override ;
    virtual void giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override ;
    virtual void giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override ;
    virtual void giveMicromorphicMatrix_dSdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override ;
    virtual void giveMicromorphicMatrix_dMdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override ;
    virtual void giveMicromorphicMatrix_dMdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override ;
    virtual void giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override ;

    virtual void giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &S, GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep) override ;
    virtual void giveFiniteStrainGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &displacementGradient, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep){;}

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

protected:
  virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MicromorphicMaterialStatus(1, domain, gp); }
                                                                            
};


} // end namespace oofem
#endif
