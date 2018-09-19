/*
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

#ifndef trabbone3d_h
#define trabbone3d_h

#include "../sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "matconst.h"
#include "matstatus.h"
#include "../sm/Materials/structuralms.h"
#include "cltypes.h"
#include "Materials/linearelasticmaterial.h"

///@name Input fields for TrabBone3D
//@{
#define _IFT_EllipticPlasticMaterial_Name "ellipticplasticmat"
#define _IFT_EllipticPlasticMaterial_sig0 "sig0"
#define _IFT_EllipticPlasticMaterial_h "h"
#define _IFT_EllipticPlasticMaterial_c "c"
#define _IFT_EllipticPlasticMaterial_f "f"


//@}

namespace oofem {
/**
 * This class implements associated Material Status to TrabBone3D (trabecular bone material).
 * It is attribute of matStatusDictionary at every GaussPoint, for which this material
 * is active.
 */
class EllipticPlasticMaterialStatus : public StructuralMaterialStatus
{
protected:
    double kappa, tempKappa, dam;
    FloatArray tempPlasticStrain, plasticStrain;



public:
    EllipticPlasticMaterialStatus(int n, Domain *d, GaussPoint *g);

    virtual ~EllipticPlasticMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    void setTempKappa(double al) { tempKappa = al; }
    void setKappa(double values) { kappa = values; }
    void setTempPlasticStrain(FloatArray &epsip) { tempPlasticStrain = epsip; }



    double giveKappa(){return kappa;}
    double giveTempKappa(){ return tempKappa;}

    const FloatArray &givePlasticStrain() const{ return plasticStrain;}
    const FloatArray &giveTempPlasDef() const{return tempPlasticStrain;}
    
    
    virtual const char *giveClassName() const { return "EllipticPlasticMaterial"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/////////////////////////////////////////////////////////////////
////////////////TRABECULAR BONE 3D///////////////////////////////
/////////////////////////////////////////////////////////////////


class EllipticPlasticMaterial : public StructuralMaterial
{
protected:
  double C, F, sig0, H;
  LinearElasticMaterial *linearElasticMaterial;
  
public:
    EllipticPlasticMaterial(int n, Domain *d);
    ~EllipticPlasticMaterial();

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }
    virtual double evaluateCurrentYieldStress(GaussPoint *gp);
    double evaluateCurrentPlasticModulus(const double kappa);

    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep);

    void constructYieldFunctionTensor(FloatMatrix &answer);


    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseMode, GaussPoint * gp,
                                               TimeStep * tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *,
                                         const FloatArray &, TimeStep *);

    virtual const char *giveInputRecordName() const { return _IFT_EllipticPlasticMaterial_Name; }
    virtual const char *giveClassName() const { return "EllipticPlasticMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

};
} //end namespace oofem
#endif
