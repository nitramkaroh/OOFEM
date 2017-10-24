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

#ifndef misesmatlogestrain_h
#define misesmatlogstrain_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/misesmat.h"
#include "Materials/linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for MisesMat
//@{
#define _IFT_MisesMatLogStrain_Name "misesmatlogstrain"
#define _IFT_MisesMatLogStrain_K "k"
#define _IFT_MisesMatLogStrain_mu "mu"
#define _IFT_MisesMatLogStrain_sig0 "sig0"
#define _IFT_MisesMatLogStrain_h "h"
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
 class MisesMatLogStrain : public MisesMat
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
    MisesMatLogStrain(int n, Domain * d);
    virtual ~MisesMatLogStrain();




    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }

    virtual const char *giveInputRecordName() const { return _IFT_MisesMatLogStrain_Name; }
    virtual const char *giveClassName() const { return "MisesMatLogStrain"; }



    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

   
    virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);
    virtual void giveCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);


virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
						GaussPoint *gp, TimeStep *tStep);
    
    virtual void giveSpatial3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep);



protected:
    


    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    void compute_dElog_dE(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &lam, const FloatMatrix &N);
    void compute_dB_dF(FloatMatrix &answer, const FloatMatrix &Btr);
    void convert_LogStiffness_2_SpatialStiffness(FloatMatrix &answer, const FloatMatrix &stiffness,MatResponseMode  mode,  GaussPoint *gp, TimeStep *tStep);
    void  give_Sigma_Delta_Product( FloatMatrix &answer, const FloatArray &vCauchy);
    void giveDLBProduct(FloatMatrix &answer, const FloatMatrix &DD, const FloatMatrix &LL, const FloatMatrix &BB);



};


class MisesMatLogStrainStatus : public MisesMatStatus
{
protected:
 



public:
    MisesMatLogStrainStatus(int n, Domain * d, GaussPoint * g);
    virtual ~MisesMatLogStrainStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual const char *giveClassName() const { return "MisesMatLogStrainStatus"; }
};
} // end namespace oofem
#endif // misesmat_h
