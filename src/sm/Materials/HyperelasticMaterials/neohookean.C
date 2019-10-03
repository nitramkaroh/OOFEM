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

#include "neohookean.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_Material(NeoHookeanMaterial);

NeoHookeanMaterial :: NeoHookeanMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{ }



void 
NeoHookeanMaterial :: giveCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
{
  StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );


  double J;
  FloatArray vFn, tempvF;
  FloatMatrix F, Fn, f, b, I, devC, p, Cauchy;
  vFn = status->giveFVector();
  Fn.beMatrixForm(vFn);
  f.beMatrixForm(vF);
  F.beProductOf(f,Fn);
  tempvF.beVectorForm(F);
  J = F.giveDeterminant();
  b.beProductTOf(F,F);
  I.resize(3,3);
  I.beUnitMatrix();
  p = I;
  p.times((lambda/J)*log(J));
  devC = b;
  devC.subtract(I);
  devC.times(mu/J);

  Cauchy = p;
  Cauchy.add(devC);

  answer.beSymVectorForm(Cauchy);

  
 // update gp
  status->letTempFVectorBe(tempvF);
  status->letTempCVectorBe(answer);


}


void 
NeoHookeanMaterial :: give3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J, l, m;
    FloatArray vF = status->giveTempFVector();
    FloatMatrix F, I1, I2;

    F.beMatrixForm(vF);
    J = F.giveDeterminant();

    // delta(i,j)*delta(k,l)
    I1.resize(6,6);
    I1.at(1,1) = I1.at(1,2) = I1.at(1,3) = 1;
    I1.at(2,1) = I1.at(2,2) = I1.at(2,3) = 1;
    I1.at(3,1) = I1.at(3,2) = I1.at(3,3) = 1;
    // delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k)     
    I2.resize(6,6);
    I2.at(1,1) = I2.at(2,2) = I2.at(3,3) = 2;
    I2.at(4,4) = I2.at(5,5) = I2.at(6,6) = 1;

    l = lambda/J;
    m = (mu - lambda*log(J))/J;

    I1.times(l);
    I2.times(m);
    answer = I1;
    answer.add(I2);
}




MaterialStatus *
NeoHookeanMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(1, this->giveDomain(), gp);
}


IRResultType
NeoHookeanMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro


    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    IR_GIVE_FIELD(ir, lambda, _IFT_NeoHookeanMaterial_lambda);
    IR_GIVE_FIELD(ir, mu, _IFT_NeoHookeanMaterial_mu);


    return IRRT_OK;
}

} // end namespace oofem
