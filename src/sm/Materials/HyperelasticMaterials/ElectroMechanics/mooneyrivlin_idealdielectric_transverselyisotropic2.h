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

#ifndef mooneyrivlin_idealdielectric_transverselyisotropic2_h
#define mooneyrivlin_idealdielectric_transverselyisotropic2_h

//#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/HyperelasticMaterials/ElectroMechanics/mooneyrivlin_idealdielectric.h"
#include "../sm/Materials/ElectroMechanics/electromechanicalms.h"
#include "../sm/Materials/HyperelasticMaterials/mooneyrivlin.h"
#include "../sm/Materials/ElectroMechanics/electromechanicalmaterialextensioninterface.h"

#include "cltypes.h"
#include "material.h"

///@name Input fields for MicromorphLEmat
//@{
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_Name "mooneyrivlin_idealdielectric_transverselyisotropic_mat2"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon1 "epsilon1"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon2 "epsilon2"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon3 "epsilon3"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon4 "epsilon4"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon5 "epsilon5"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_mu1 "mu1"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_mu2 "mu2"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_mu3 "mu3"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_lambda "lambda"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_n "n"

//@}

namespace oofem {
/**
 * Electromechanical coupling considering Ideal dielectric Mooney-Rivlin material
 */
  class MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2 :  public Material, public ElectroMechanicalMaterialExtensionInterface_3Fields
{
protected:
  MooneyRivlinMaterial *hyperelasticMaterial;
  double lambda, mu1, mu2, mu3;
  double iEps1, iEps2, iEps3, iEps4, iEps5;
  FloatArray N;
  double alpha, beta;
  
public:
   MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2(int n, Domain * d);
    virtual ~MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2(){;}

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_Name; }
    virtual const char *giveClassName() const { return "MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    
    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == ElectroMechanicalMaterialExtensionInterface_3FieldsType) {
            return static_cast< ElectroMechanicalMaterialExtensionInterface_3Fields* >(this);
        } else {
            return NULL;
        }
    }
    
    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) {    OOFEM_ERROR("Shouldn't be called."); }

  
    virtual void give_FirstPKStressVector_ElectricalFieldVector_3d(FloatArray &P, FloatArray &E, GaussPoint *gp, const FloatArray &F, const FloatArray &D, TimeStep *tStep) override;

    
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
    int computeAcousticTensorMinEigenvalue(GaussPoint *gp, TimeStep *tStep);


    /// Invariant FN \cdot FN
    double compute_I4(const FloatMatrix &F);
    /// Invariant HN \cdot HN
    double compute_I5(const FloatMatrix &F);
    /// Invariant D_0 \cdot  D_0
    double compute_I6(const FloatArray &D);
    /// Invariant FD \cdot FD
    double compute_I7(const FloatMatrix &F, const FloatArray &D);
    /// Invariant HD \cdot HD
    double compute_I8(const FloatMatrix &F, const FloatArray &D);
    /// (N\cdotD)^2
    double compute_I9(const FloatArray &D);
    /// (FN \cdot FD)^2
    double compute_I10(const FloatMatrix &F, const FloatArray &D_0);


    void compute_dI1dF(FloatArray &answer, const FloatMatrix &F);
    void compute_dI2dF(FloatArray &answer, const FloatMatrix &F);
    
    void compute_dI4dF(FloatArray &answer, const FloatMatrix &F);
    void compute_dI5dF(FloatArray &answer, const FloatMatrix &F);
    // dI6dF is zero
    void compute_dI7dF(FloatArray &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_dI8dF(FloatArray &answer, const FloatMatrix &F, const FloatArray &D);
    // dI9dF is zero
    void compute_dI10dF(FloatArray &answer, const FloatMatrix &F, const FloatArray &D);
    ///////////////////////////////////////////////////////////////
    
    // dI4dF and dI5dF are zero
    void compute_dI6dD(FloatArray &answer, const FloatArray &D);
    void compute_dI7dD(FloatArray &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_dI8dD(FloatArray &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_dI9dD(FloatArray &answer, const FloatArray &D);
    void compute_dI10dD(FloatArray &answer, const FloatMatrix &F, const FloatArray &D);
    ///////////////////////////////////////////////////////////////
    void compute_d2I1dF2(FloatMatrix &answer, const FloatMatrix &F);
    void compute_d2I2dF2(FloatMatrix &answer, const FloatMatrix &F);
    void compute_d2I4dF2(FloatMatrix &answer, const FloatMatrix &F);
    void compute_d2I5dF2(FloatMatrix &answer, const FloatMatrix &F);
    // d2I6dF2 is zero
    void compute_d2I7dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I8dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I9dF2(FloatMatrix &answer, const FloatArray &D);
    void compute_d2I10dF2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    ////////////////////////////////////////////////////////////////
    void compute_d2I6dD2(FloatMatrix &answer, const FloatArray &D);
    void compute_d2I7dD2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I8dD2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I9dD2(FloatMatrix &answer, const FloatArray &D);
    void compute_d2I10dD2(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    ////////////////////////////////////////////////////////////////
    void compute_d2I7dFdD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I8dFdD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_d2I10dFdD(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);    

    void compute_dI4dF_num(FloatArray &answer, const FloatMatrix &F);
    void compute_dI7dF_num(FloatArray &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_dI10dF_num(FloatArray &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_dI7dD_num(FloatArray &answer, const FloatMatrix &F, const FloatArray &D);
    void compute_dI10dD_num(FloatArray &answer, const FloatMatrix &F, const FloatArray  &D);
    void compute_d2I4dF2_num(FloatMatrix &answer,const FloatMatrix &F);
    void compute_d2I7dF2_num(FloatMatrix &answer,const FloatMatrix &F, const FloatArray &D);
    void compute_d2I10dF2_num(FloatMatrix &answer,const FloatMatrix &F, const FloatArray &D);
    void  compute_d2I7dFdD_num(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void  compute_d2I10dFdD_num(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);
    void  compute_d2I7dD2_num(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &D);


    void compute_dPdF_dEdF(FloatMatrix &dPdF,FloatMatrix &dEdF, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep);
    void compute_dPdF_dEdF_fake(FloatMatrix &dPdF,FloatMatrix &dEdF, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep);
    void compute_dPdD_dEdD(FloatMatrix &dPdF,FloatMatrix &dEdF, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep);


    void computePseudoInverse(FloatMatrix &iDdd, FloatMatrix &b);

    
};
 

} // end namespace oofem
#endif
