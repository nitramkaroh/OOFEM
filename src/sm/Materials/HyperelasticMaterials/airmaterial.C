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

#include "airmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "simopistermat.h"
#include "mooneyrivlin.h"
#include "domain.h"
#include "function.h"


namespace oofem {
  REGISTER_Material(AirMaterial);

AirMaterial :: AirMaterial(int n, Domain *d) :StructuralMaterial(n, d)
{
     
}


AirMaterial :: ~AirMaterial()
{
  delete regularizationMaterial;
}

  

void
AirMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{
    FloatArray vH;
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );    
    //store deformation gradient into matrix
    this->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    double factor = p * this->giveDomain()->giveFunction(ltf)->evaluate(tStep, VM_Total);
    //1.PK
    answer.zero();    
    answer.add( factor, vH);
    // regularization part
    FloatArray vP;
    regularizationMaterial->giveFirstPKStressVector_3d(vP, gp, vF, tStep);
    answer.add(eps, vP);
    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
}

void
AirMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    answer.resize(9,9);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    //deformation gradient from the status
    FloatArray vF;
    vF = status->giveTempFVector();
    //
    FloatMatrix Fx;
    this->compute_tensor_cross_product_tensor(Fx, vF);
    double factor = p * this->giveDomain()->giveFunction(ltf)->evaluate(tStep, VM_Total);
    //
    FloatMatrix D;
    regularizationMaterial->give3dMaterialStiffnessMatrix_dPdF(D, mode, gp, tStep);
    //
    answer.zero();
    answer.add(factor, Fx);
    answer.add(eps, D); 
}



  


MaterialStatus *
AirMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


IRResultType
AirMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    //
    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, p, _IFT_AirMaterial_p);
    IR_GIVE_FIELD(ir, ltf, _IFT_AirMaterial_ltf);
    IR_GIVE_FIELD(ir, eps, _IFT_AirMaterial_eps);

    int regularizationMaterialType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, regularizationMaterialType, _IFT_AirMaterial_rmt);

    if ( regularizationMaterialType == 0 ) {
      regularizationMaterial = new SimoPisterMaterial(this->giveNumber(), this->giveDomain());
    } else if ( regularizationMaterialType == 1 ) {
      regularizationMaterial = new MooneyRivlinMaterial(this->giveNumber(), this->giveDomain());
    } else {
      OOFEM_WARNING("Unknown hyperelasticmaterial type %d", regularizationMaterialType);
      return IRRT_BAD_FORMAT;
    }
    //
    regularizationMaterial->initializeFrom(ir);
    if ( result != IRRT_OK ) {
      return result;
    }
    //   
    return IRRT_OK;
}



} // end namespace oofem
