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

#include "../sm/Materials/Micromorphic/GradientPolyconvex/stvenantkirchhoff_incompressible.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"
#include "engngm.h"

namespace oofem {
  REGISTER_Material(StVenantKirchhoffIncompressibleMaterial);

StVenantKirchhoffIncompressibleMaterial :: StVenantKirchhoffIncompressibleMaterial(int n, Domain *d) :IsotropicLinearElasticMaterial(n, d)
{
}

StVenantKirchhoffIncompressibleMaterial :: ~StVenantKirchhoffIncompressibleMaterial()
{ }



void
StVenantKirchhoffIncompressibleMaterial :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}

void
StVenantKirchhoffIncompressibleMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{
  StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    FloatMatrix F;
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    // compute jacobian
    double J = F.giveDeterminant();
    // St. Venant-Kirchhoff material model
    IsotropicLinearElasticMaterial :: giveFirstPKStressVector_3d(answer, gp, vF, tStep);

    FloatArray vPm(9);
    vPm.zero();
    double lnJ = log(J);
    this->compute_dJ_dF(vPm, F);

    if(tStep->giveIntrinsicTime() > 280) { 
      vPm.times(K * lnJ / J);
      answer.add(vPm);
    }


    /*    if(J <= 0.) {
      this->giveDomain()->giveEngngModel()->setAnalysisCrash(true);
      }*/

 
    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
}



void
StVenantKirchhoffIncompressibleMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J;
    FloatArray vF;
    FloatMatrix F, dSdE, invF, invFt;
    FloatMatrix dInvF_dF, iFtxiFt;
    this->give3dMaterialStiffnessMatrix(dSdE, mode, gp, tStep);
    this->give_dPdF_from(dSdE, answer, gp, _3dMat);

    //deformation gradient from the status
    vF = status->giveTempFVector();
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    J = F.giveDeterminant();

    invF.beInverseOf(F);
    invFt.beTranspositionOf(invF);

    double lnJ = log(J);
    
    this->compute_dInvFt_dF(dInvF_dF, invF);
    this->compute_dyadic_product(iFtxiFt, invFt, invFt);
    
    iFtxiFt.times(K);
    dInvF_dF.times(K * lnJ);
    if(tStep->giveIntrinsicTime() > 280) { 
      answer.add(iFtxiFt);
      answer.add(dInvF_dF); 
    }


}
  


  
  

IRResultType
StVenantKirchhoffIncompressibleMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    IsotropicLinearElasticMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, K, _IFT_StVenantKirchhoffIncompressibleMaterial_K);
    
    
    return IRRT_OK;
}



} // end namespace oofem
