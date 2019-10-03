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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "ogdenmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_Material(OgdenNematicElastomerMaterial);

OgdenNematicElastomerMaterial :: OgdenNematicElastomerMaterial(int n, Domain *d) : OgdenMaterial(n, d)
{ }


 

void
OgdenNematicElastomerMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{

    FloatArray eVals;
    FloatMatrix P, F, S, C, eVecs;
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    //
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    if(!qcEnvelop) {    
      this->giveSecondPKStressVector_3d(S, C, eVals, eVecs, tStep);
      // add the terms according to nematic elastomer
    } else {
      if(lambda.at(3) <= pow(a,1./6.)) {
	S.resize(3,3);
	S.zero();
      } else if(1./lambda.at(3) >= power(a, -0.5) * lambda.at(3)) {
	this->giveSecondPKStressVector_3d(S, C, eVals, eVecs, tStep);
      } else {
	this->giveSecondPKStressVectorModified_3d(S, C, eVals, eVecs, tStep);
      }      
    }
    
    P.beProductOf(F,S);
    answer.beVectorForm(P);
    // update gp
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
}



void
OgdenNematicElastomerMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J, lnJ;
    FloatArray vF, eVals;
    FloatMatrix P, F, S, C, invC, iCiC, dSdE;
    FloatMatrix eVecs, invF, invFt, d2I1dF2, d2I2dF2, dinvF_dF, iFtxiFt, dInvF_dF;
    
    //deformation gradient from the status
    vF = status->giveTempFVector();
    //store deformation gradient into matrix
    F.beMatrixForm(vF);
    J = F.giveDeterminant();
    lnJ = log(J);
    //
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);

    if(!qcEnvelop) {    
      this->give3dMaterialStiffnessMatrix_dSdE(dSdE, mode, gp, tStep);
      // add the terms according to nematic elastomer
    } else {
      if(lambda.at(3) <= pow(a,1./6.)) {
	S.resize(3,3);
	S.zero();
      } else if(1./lambda.at(3) >= power(a, -0.5) * lambda.at(3)) {
	this->give3dMaterialStiffnessMatrix_dSdE(dSdE, mode, gp, tStep);
      } else {
	this->give3dMaterialStiffnessMatrixModified_dSdE(dSdE, mode, gp, tStep);
      }      
    }
   
    this->give_dPdF_from(dSdE, answer, gp, _3dMat);
    
}



  


MaterialStatus *
OgdenNematicElastomerMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


IRResultType
OgdenNematicElastomerMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    result = OgdenMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, a, _IFT_OgdenNematicElastomerMaterial_a);
    this->qcEnvelop = ir->hasField(_IFT_OgdenNematicElastomerMaterial_qce);
    
    return IRRT_OK;
}

} // end namespace oofem
