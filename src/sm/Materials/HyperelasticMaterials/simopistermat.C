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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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

#include "simopistermat.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(SimoPisterMaterial);

SimoPisterMaterial::SimoPisterMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }


void
SimoPisterMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{
    double J;
    FloatMatrix F;

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    // compute jacobian and its logarith
    J = F.giveDeterminant();
    // Derivatives of the invariants
    FloatArray dI1_dF, vH;
    this->compute_dI1_C_dF(dI1_dF, F);
    //
    this->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    //1.PK
    answer.zero();
    answer.add(0.5 * mu, dI1_dF);
    answer.add(-mu / J, vH);
    answer.add(lambda * ( J - 1. ), vH);
    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
}

void
SimoPisterMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J;
    FloatArray vF, vH;
    FloatMatrix F, d2I1dF2, d2JdF2;
    //deformation gradient from the status
    vF = status->giveTempFVector();
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    J = F.giveDeterminant();
    //
    this->compute_2order_tensor_cross_product( vH, vF, vF );
    vH.times(0.5);
    //
    FloatMatrix HxH;
    HxH.beProductTOf(vH,vH);     
    //
    FloatMatrix Fx;
    this->compute_tensor_cross_product_tensor(Fx, vF);
    //
    this->compute_d2I1_C_dF2(d2I1dF2, F);
    //
    answer.add(0.5 * mu, d2I1dF2);
    answer.add(mu / J / J + lambda, HxH);
    answer.add(lambda * (J-1)  - mu / J, Fx);
}



  


MaterialStatus *
SimoPisterMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


IRResultType
SimoPisterMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    //
    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, mu, _IFT_SimoPisterMaterial_mu);
    IR_GIVE_FIELD(ir, lambda, _IFT_SimoPisterMaterial_lambda);

    return IRRT_OK;
}

} // end namespace oofem
