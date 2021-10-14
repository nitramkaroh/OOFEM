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
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon1 "eps1"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon2 "eps2"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon3 "eps3"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon4 "eps4"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_epsilon5 "eps5"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_mu1 "mu1"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_mu2 "mu2"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_mu3 "mu3"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_lambda "lambda"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_n "n"

#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_alpha "alpha"
#define _IFT_MooneyRivlin_IdealDielectric_TransverselyIsotropicMaterial2_beta "beta"

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
  double a_I8pol = 1., b_I8pol = 1., c_I8pol = 1.;
  double a_K2_Cinf_pol = 1., b_K2_Cinf_pol = 1.;
  double a_K2_Dinf_pol = 1., b_K2_Dinf_pol = 1.;
  
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
    ///////////////////////////////////////
    void compute_dPdD_dEdD(FloatMatrix &dPdD,FloatMatrix &dEdD, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep);
    void compute_dPdF_dEdF(FloatMatrix &dPdF,FloatMatrix &dEdF, const FloatMatrix &F, const FloatArray &D, GaussPoint *gp, TimeStep *tStep);




    /// Invariant F:F
    double compute_I1(double J, const FloatArray &vH, const FloatArray &vF);
    /// Invariant H:H
    double compute_I2(double J, const FloatArray &vH, const FloatArray &vF);
    /// Invariant J^{-2/3}F:F
    double compute_I1_dev(double J, const FloatArray &vH, const FloatArray &vF);
    /// Invariant J^{-4/3}F:F
    double compute_I2_dev(double J, const FloatArray &vH, const FloatArray &vF);
    /// Invariant (I2_dev)^{3/2}
    double compute_I2_dev_pol(double J, const FloatArray &vH, const FloatArray &vF);
    /// Invariant det F = 1/3 F:H
    double compute_I3(double J, const FloatArray &vH, const FloatArray &vF);
    /// Invariant FN \cdot FN
    double compute_I4(double J, const FloatArray &vH, const FloatArray &vF);
    /// Invariant HN \cdot HN
    double compute_I5(double J, const FloatArray &vH, const FloatArray &vF);
    /// Invariant D_0 \cdot  D_0
    double compute_I6(const FloatArray &D);
    /// Invariant FD \cdot FD
    double compute_I7(double J, const FloatArray &vH, const FloatArray &vF,   const FloatArray &D);
    /// Invariant HD \cdot HD
    double compute_I8(double J, const FloatArray &vH, const FloatArray &vF,   const FloatArray &D);
    double compute_I8_pol(double J, const FloatArray &vH, const FloatArray &vF,   const FloatArray &D);
    // D\cdot N
    double compute_K1_Cinf(const FloatArray &D);
    // d \cdot FN
    double compute_K2_Cinf(double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    double compute_K2_Cinf_pol(double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    double compute_K1_Dinf(const FloatArray &D);
    double compute_K2_Dinf(double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    double compute_K2_Dinf_pol(double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    //derivatives
    void compute_dI1_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_dI2_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF);  
    void compute_dI3_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_dI1dev_dF(FloatArray &dIdF, const double J, const FloatArray &H, const FloatArray &F);
    void compute_dI2dev_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF);   
    void compute_dI2devpol_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_dI4_dF(FloatArray &dIdF, const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_dI5_dF(FloatArray &dIdF, const double J, const FloatArray &H, const FloatArray &F);
    void compute_dI6_dD(FloatArray &dIdF, const FloatArray &D);
    void compute_dI7_dF_dD(FloatArray &dIdF, FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_dI8_dF_dD(FloatArray &dIdF, FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_dI8pol_dF_dD(FloatArray &dIdF, FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_dK1Cinf_dD(FloatArray &answer, const FloatArray &D);
    void compute_dK2Cinf_dF_dD(FloatArray &dIdF,FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_dK2Cinfpol_dF_dD(FloatArray &dIdF,FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_dK1Dinf_dD(FloatArray &answer, const FloatArray &D);
    void compute_dK2Dinf_dF_dD(FloatArray &dIdF,FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_dK2Dinfpol_dF_dD(FloatArray &dIdF,FloatArray &dIdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void compute_d2I1_dF2(FloatMatrix &d2IdF2,const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_d2I2_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF);   
    void compute_d2I3_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &H, const FloatArray &F);
    void compute_d2I1dev_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_d2I2dev_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_d2I2devpol_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_d2I4_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_d2I5_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF);
    void compute_d2I6_dD2(FloatMatrix &answer, const FloatArray &D);
    void compute_d2I7_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2I7_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2I7_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2I8_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2I8_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2I8_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2I8pol_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2I8pol_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2I8pol_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K1Dinf_dD2(FloatMatrix &d2IdD2, const FloatArray &D);
    void compute_d2K2Cinf_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Cinf_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Cinf_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Cinfpol_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Cinfpol_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Cinfpol_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Dinf_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Dinf_dD2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Dinf_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Dinfpol_dF2(FloatMatrix &d2IdF2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Dinfpol_dD2(FloatMatrix &d2IdD2, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    void compute_d2K2Dinfpol_dFdD(FloatMatrix &d2IdFdD, const double J, const FloatArray &vH, const FloatArray &vF,  const FloatArray &D);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_dI1dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_dI2dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_dI1devdF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_dI2devdF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_dI2devpoldF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_dI4dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_dI5dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_dI7dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dI8dF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dI8poldF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dK2CdF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dK2CpoldF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dK2DdF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dK2DpoldF_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_dI6dD_num(FloatArray &answer,  const FloatArray &D);
void compute_dI7dD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dI8dD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dI8poldD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dK1CdD_num(FloatArray &answer, const FloatArray &D);
void compute_dK2CdD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dK2CpoldD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_dK1DdD_num(FloatArray &answer, const FloatArray &D);
void compute_dK2DdD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D); 
void compute_dK2DpoldD_num(FloatArray &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_d2I1dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_d2I1devdF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF);
 void compute_d2I2dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_d2I2devdF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_d2I2devpoldF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_d2I4dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_d2I5dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF);
void compute_d2I7dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2I8dF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2I8poldF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2CdF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2CpoldF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2DdF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2DpoldF2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_d2I6dD2_num(FloatMatrix &answer, const FloatArray &D);
void compute_d2I7dD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2I8dD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2I8poldD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2CdD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2CpoldD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K1DdD2_num(FloatMatrix &answer,  const FloatArray &D);
void compute_d2K2DdD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2DpoldD2_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
///////////////////////////////////////////////////////////////////////////////////////////////////
void compute_d2I7dFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2I8dFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2I8poldFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2CdFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
 void compute_d2K2CpoldFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2DdFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);
void compute_d2K2DpoldFdD_num(FloatMatrix &answer, const double J, const FloatArray &vH, const FloatArray &vF, const FloatArray &D);

    
};
 

} // end namespace oofem
#endif
