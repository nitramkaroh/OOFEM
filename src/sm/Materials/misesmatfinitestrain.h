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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#ifndef misesmatfinitestrain_h
#define misesmatfinitestrain_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/incompressiblematerialextensioninterface.h"
#include "Materials/linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for MisesMat
//@{
#define _IFT_MisesMatFiniteStrain_Name "misesmatfinitestrain"
#define _IFT_MisesMatFiniteStrain_K "k"
#define _IFT_MisesMatFiniteStrain_mu "mu"
#define _IFT_MisesMatFiniteStrain_sig0 "sig0"
#define _IFT_MisesMatFiniteStrain_h "h"
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * This class implements an isotropic elastoplastic material
 * with Mises yield condition, associated flow rule
 * and linear isotropic hardening.
 *
 * It differs from other similar materials (such as J2Mat)
 * by implementation - here we use the radial return, which
 * is the most efficient algorithm for this specific model.
 * Also, an extension to large strain will be available.
 * 
 * The model also exemplifies how to implement non-3d material modes, in this case 1D, 
 * by overloading the default implementations that iterates over the 3D method.
 */
 class MisesMatFiniteStrain : public StructuralMaterial, public IncompressibleMaterialExtensionInterface
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial;

    /// Elastic shear modulus.
    double G;

    /// Elastic bulk modulus.
    double K;

    /// Hardening modulus.
    double H;

    /// Initial (uniaxial) yield stress.
    double sig0;

    double mu;
public:
    MisesMatFiniteStrain(int n, Domain * d);
    virtual ~MisesMatFiniteStrain();


    virtual void computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }

    virtual const char *giveInputRecordName() const { return _IFT_MisesMatFiniteStrain_Name; }
    virtual const char *giveClassName() const { return "MisesMatFiniteStrain"; }


    virtual Interface *giveInterface(InterfaceType t) {
        if ( t ==  IncompressibleMaterialExtensionInterfaceType ) {
            return static_cast<  IncompressibleMaterialExtensionInterface * >(this);
        } else {
            return NULL;
        }
    }


    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep){OOFEM_ERROR("This material can be used only in the large strain mode");}
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep){OOFEM_ERROR("This material can be used only in the large strain mode");} 
  
  

    virtual void giveCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);

    virtual void give3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep);



    double computePressure(double J);
    virtual void giveDeviatoricCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);
    virtual void giveVolumetricCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const double J);
    virtual void giveInitialStiffnessMatrix_Cauchy(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    virtual void givePressure3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, double J);
    virtual void giveVolumetric3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, double J);
    virtual void giveDeviatoric3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);



protected:
    

   double  giveSl(double lambda_alpha,double lambda_beta, double sigma_alpha, double sigma_beta,double J); 

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};


class MisesMatFiniteStrainStatus : public StructuralMaterialStatus
{
protected:
 
   /// inversion of right plastic Cauchy-Green tensor (initial).
    FloatArray invCp;

    /// inversion of right plastic Cauchy-Green tensor (final).
    FloatArray tempInvCp;


    /// Cumulative plastic strain (initial).
    double kappa;

    /// Cumulative plastic strain (final).
    double tempKappa;

    /// deviatoric part of Cauchy stress (initial)
    FloatArray devCauchyStressVector;
    /// deviatoric part of Cauchy stress (final)
    FloatArray tempDevCauchyStressVector;

    /// deviatoric part of Cauchy stress (initial)
    FloatArray volCauchyStressVector;
    /// deviatoric part of Cauchy stress (final)
    FloatArray tempVolCauchyStressVector;


    /// deviatoric trial stress - needed for tangent stiffness.
    FloatArray trialStressDev;
    FloatArray trialLambda;



public:
    MisesMatFiniteStrainStatus(int n, Domain * d, GaussPoint * g);
    virtual ~MisesMatFiniteStrainStatus();

    double giveCumulativePlasticStrain() { return kappa; }
    double giveTempCumulativePlasticStrain() { return tempKappa; }
    const FloatArray & giveInvLeftCauchyGreenPl() { return invCp; }
    const FloatArray & giveTempInvLeftCauchyGreenPl() { return tempInvCp; }
    const FloatArray &giveTrialStressDev() { return trialStressDev; }    
    const FloatArray &giveTrialLambda() { return trialLambda; }    
    const FloatArray &giveTempDeviatoricCauchyStressVector(){return tempDevCauchyStressVector;}
    const FloatArray &giveTempVolumetricCauchyStressVector(){return tempVolCauchyStressVector;}


    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }
    void letTempInvLeftCauchyGreenPlBe(const FloatArray &values) { tempInvCp = values; }
    void letTempDeviatoricCauchyStressVectorBe(const FloatArray &values){tempDevCauchyStressVector = values;}
    void letTempVolumetricCauchyStressVectorBe(const FloatArray &values){tempVolCauchyStressVector = values;}
    void letTrialStressDevBe(const FloatArray &values) { trialStressDev = values; }
    void letTrialLambdaBe(const FloatArray &values) { trialLambda = values; }


    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "MisesMatFiniteStrainStatus"; }
};
} // end namespace oofem
#endif // misesmat_h
