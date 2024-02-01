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

#ifndef airmaterialcontact_h
#define airmaterialcontact_h


#include "../sm/Materials/Micromorphic/secondgradientmaterialextensioninterface.h"
#include "../sm/Materials/Micromorphic/secondgradientms.h"
#include "../sm/Materials/HyperelasticMaterials/airmaterial.h"
#include "cltypes.h"


///@name Input fields for MicromorphLEmat
//@{
#define _IFT_AirMaterialContactRot_Name "airmaterialcontactrot"
#define _IFT_AirMaterialContactRot_a "a"
//@}

namespace oofem {

  class AirMaterialContactRotStatus : public SecondGradientMaterialStatus
 {
 public:
   AirMaterialContactRotStatus(int n, Domain *d, GaussPoint *g);  
   ~AirMaterialContactRotStatus(){;}
 };

/**
 * MicromorphicLinearElasticMaterial
 */
class AirMaterialContactRot : public AirMaterial, SecondGradientMaterialExtensionInterface
{
protected:
  double a;
public:
    AirMaterialContactRot(int n, Domain * d);
    virtual ~AirMaterialContactRot();

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_AirMaterialContactRot_Name; }
    virtual const char *giveClassName() const { return "AirMaterialContactRot"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == SecondGradientMaterialExtensionInterfaceType ) {
            return static_cast< SecondGradientMaterialExtensionInterface * >(this);
        } else {
            return NULL;
        }
    }

    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual void giveSecondGradientMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveSecondGradientMatrix_dPdG(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveSecondGradientMatrix_dMdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveSecondGradientMatrix_dMdG(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    /// micromorhpic stresses
    virtual void giveGeneralizedStressVectors (FloatArray &vP, FloatArray &vM, const FloatArray &vF, const FloatArray &vG, GaussPoint *gp, TimeStep *tStep)
  {
    this->giveFiniteStrainGeneralizedStressVectors(vP, vM, vF, vG, gp, tStep);
  }
  virtual void giveFiniteStrainGeneralizedStressVectors(FloatArray &vP, FloatArray &vM, const FloatArray &vF, const FloatArray &vG, GaussPoint *gp, TimeStep *tStep);
  virtual void giveFiniteStrainGeneralizedStressVectors_3d (FloatArray &vP, FloatArray &vM, const FloatArray &vF, const FloatArray &vG, GaussPoint *gp, TimeStep *tStep);
  virtual void giveFiniteStrainGeneralizedStressVectors_PlaneStrain (FloatArray &vP, FloatArray &vM, const FloatArray &vF, const FloatArray &vG, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new AirMaterialContactRotStatus(1, domain, gp); }

                                                                     
};


 
 

} // end namespace oofem
#endif
