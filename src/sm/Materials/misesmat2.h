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

#ifndef misesmat_h
#define misesmat_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "stressvector.h"

///@name Input fields for MisesMat
//@{
#define _IFT_MisesMat2_Name "misesmat2"
#define _IFT_MisesMat2_sig0 "sig0"
#define _IFT_MisesMat2_hi "hi"
#define _IFT_MisesMat2_hk "hk"
#define _IFT_MisesMat2_y "y"
#define _IFT_MisesMat2_expD "expd"
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
 */
class MisesMat2 : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial;

    /// Elastic shear modulus.
    double G;

    /// Elastic bulk modulus.
    double K;

   enum HardeningType { EST_Linear, EST_NonLinear };



    /// Initial (uniaxial) yield stress.
    double sig0;

	
    /// Isotropic Hardening modulus.
    double Hi;
	double Y;
	double expD;

	/// Kinematic Hardening modulus.
    double Hk;



public:
    MisesMat2(int n, Domain *d);
    virtual ~MisesMat2();

    void performSSPlasticityReturn(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain);  
void performSSPlasticityReturn_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain);
	void performLSPlasticityReturn(FloatMatrix &answer, GaussPoint *gp, const FloatMatrix &Bbar);  	
	virtual void computeTrialStress(StressVector &answer, GaussPoint *gp, const StrainVector &elStrainDev);
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }
    virtual const char *giveInputRecordName() const { return _IFT_MisesMat2_Name; }
    virtual const char *giveClassName() const { return "MisesMat2"; }
    /// Returns a reference to the basic elastic material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }
	    
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

	void giveMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);


    
    void give_dPdF_from_dCde(FloatMatrix &stiffness, FloatMatrix &answer, FloatMatrix &invF, FloatMatrix &kirchhoffStress, MaterialMode matMode);


   virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mode, GaussPoint *gp,
                                          TimeStep *tStep)
   {this -> giveMaterialStiffnessMatrix(answer,mode,gp,tStep);}
   
   
	
	virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

	virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    virtual void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }

    
	void giveFirstPKStressVector(FloatArray &answer,GaussPoint *gp,const FloatArray &reducedvF,TimeStep *atTime);    
	virtual void giveFirstPKStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
	{this->giveFirstPKStressVector(answer,gp,vF,tStep);}
	virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
	{this->giveFirstPKStressVector(answer,gp,vF,tStep);}
    

	virtual void giveCauchyStressVector_3d(FloatArray &answer,GaussPoint *gp,const FloatArray &incvF,TimeStep *atTime);
	virtual void giveMaterialStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);


	double computePressure(double J);

	
	virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,MatResponseMode mode,GaussPoint *gp,TimeStep *tStep);
	
   
	
	
protected:
    //virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep); 
	double computeYieldStress(double kappa);
	double computeYieldStressPrime(double kappa);

	double computeBackStressModulus(double kappa);
	double computeBackStressModulusPrime(double kappa);


	void  givePinvID(FloatMatrix &answer, MaterialMode matMode);
	void  giveDeltaVector(FloatArray &answer, MaterialMode matMode);


    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};

//=============================================================================


class MisesMat2Status : public StructuralMaterialStatus
{
protected:
    /// Plastic strain (initial).
    FloatArray plasticStrain;

    /// Plastic strain (final).
    FloatArray tempPlasticStrain;

    /// Deviatoric trial stress - needed for tangent stiffness.
    FloatArray trialStressD;

    /// Back stress(initial)
    FloatArray backStress;

    /// Back stress(final)
    FloatArray tempBackStress;

    /// Cumulative plastic strain (initial).
    double kappa;

    /// Cumulative plastic strain (final).
    double tempKappa;


    /// Left Cauchy-Green tensor(initial - only for large strain model)
    FloatMatrix leftCauchyGreenMatrix;
    /// Left Cauchy-Green tensor(final - only for large strain model)
    FloatMatrix tempLeftCauchyGreenMatrix;
    
public:
    MisesMat2Status(int n, Domain *d, GaussPoint *g);
    virtual ~MisesMat2Status();

    void givePlasticStrain(FloatArray &answer) { answer = plasticStrain; }
    void giveTrialStressDev(FloatArray &answer) { answer = trialStressD; }
    void giveBackStress(FloatArray &answer) { answer = backStress; } 
    void giveTempBackStress(FloatArray &answer) { answer = tempBackStress; } 
    void giveLeftCauchyGreenMatrix(FloatMatrix &answer) { answer = leftCauchyGreenMatrix; }
    void giveTempLeftCauchyGreenMatrix(FloatMatrix &answer) { answer = tempLeftCauchyGreenMatrix; }
    
    double giveCumulativePlasticStrain() { return kappa; }
    double giveTempCumulativePlasticStrain() { return tempKappa; }
    
    void letTempPlasticStrainBe(FloatArray &values) { tempPlasticStrain = values; }
    void letTrialStressDevBe(FloatArray &values) { trialStressD = values; }
    void letTempBackStressBe(FloatArray &values) { backStress = values; }
    void letTempLeftCauchyGreenMatrixBe(FloatMatrix &values) { tempLeftCauchyGreenMatrix = values; }
    
    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }
    
    



    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

};
} // end namespace oofem
#endif // misesmat_h
