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

#ifndef mooneyrivlin_idealdielectric_h
#define mooneyrivlin_idealdielectric_h


#include "../sm/Materials/ElectroMechanics/electromechanicalmaterialextensioninterface.h"
#include "../sm/Materials/ElectroMechanics/electromechanicalms.h"
#include "../sm/Materials/HyperelasticMaterials/mooneyrivlin.h"
#include "cltypes.h"
#include "material.h"
#include "../sm/Materials/ElectroMechanics/electromechanicalms.h"

///@name Input fields for MicromorphLEmat
//@{
#define _IFT_MooneyRivlin_IdealDielectricMaterial_Name "mooneyrivlin_idealdielectricmat"
#define _IFT_MooneyRivlin_IdealDielectricMaterial_epsilon_m "epsilon_m"
#define _IFT_MooneyRivlin_IdealDielectricMaterial_epsilon_f "epsilon_f"
#define _IFT_MooneyRivlin_IdealDielectricMaterial_Cm "cm"

//@}

namespace oofem {
/**
 * Electromechanical coupling considering Ideal dielectric Mooney-Rivlin material
 */
  class MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial : MooneyRivlin_IdealDielectricMaterial
{
protected:
  MooneyRivlinMaterial *hyperelasticMaterial;
  double Cm;
  double epsilon_m;
  double epsilon_f;
  
public:
   MooneyRivlin_IdealDielectricMaterial(int n, Domain * d);
    virtual ~MooneyRivlin_IdealDielectricMaterial(){;}

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_MooneyRivlin_IdealDielectricMaterial_Name; }
    virtual const char *giveClassName() const { return "MooneyRivlin_IdealDielectricMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    
    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == ElectroMechanicalMaterialExtensionInterfaceType) {
            return static_cast< ElectroMechanicalMaterialExtensionInterface* >(this);
        } else {
            return NULL;
        }
    }
    
    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) {    OOFEM_ERROR("Shouldn't be called."); }

  
    virtual void give_FirstPKStressVector_ElectricalFieldVector_3d(FloatArray &P, FloatArray &E, GaussPoint *gp, const FloatArray &F, const FloatArray &D, TimeStep *tStep){;}

    
    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep);

    virtual void give3dMaterialStiffnessMatrix_dPdD(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep);
	
    virtual void give3dMaterialStiffnessMatrix_dEdD(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep);

    
    
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new ElectroMechanicalMaterialStatus(1, domain, gp); }

    double computeI4(const FloatMatrix &F);
    double computeI7(const FloatMatrix &F, const FloatArray &D);
    double computeI10(const FloatMatrix &F, const FloatArray &D_0);
    void compute_dI4dF(FloatMatrix &answer, const FloatMatrix &F);
    void compute_dI7dF(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_dI7dD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_dI10dF(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_dI10dD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I4dF2(FloatMatrix &answer, const FloatMatrix &F);
    void compute_d2I7dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I10dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I7dD2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I7dFdF(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I7dFdF(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);    
                                                                     
};


 
 

} // end namespace oofem
#endif
