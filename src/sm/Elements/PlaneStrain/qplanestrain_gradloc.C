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

#include "../sm/Elements/PlaneStrain/qplanestrain_gradloc.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "classfactory.h"
#include "fei2dquadquad.h"


namespace oofem {
REGISTER_Element(QPlaneStrainGradLoc);


QPlaneStrainGradLoc :: QPlaneStrainGradLoc(int n, Domain *aDomain) :
  QPlaneStrain(n, aDomain), BaseMicromorphicElement()
{

}

void
QPlaneStrainGradLoc :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {D_u, D_v, D_u, D_v};
}


void
QPlaneStrainGradLoc :: giveDofManDofIDMask_u(IntArray &answer)
{
  answer = {D_u, D_v};
}


void
QPlaneStrainGradLoc :: giveDofManDofIDMask_m(IntArray &answer)
{
  answer = {D_u, D_v};
}


void
QPlaneStrainGradLoc ::  postInitialize() 
{
  BaseMicromorphicElement :: postInitialize();
  QPlaneStrain :: postInitialize();
  locationArray_u = {1,2,5,6, 9,10,13,14,17,18,21,22,25,26,29,30};
  locationArray_m = {3,4,7,8,11,12,15,16,19,20,23,24,27,28,31,32};
  /*locationArray_u = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};  
  locationArray_m = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};  
  */
}



void
QPlaneStrainGradLoc :: computeMicromorphicVars(FloatArray &micromorphicVar,FloatArray &micromorphicVarGrad, IntArray IdMask_m, GaussPoint *gp, TimeStep *tStep)
{
  FloatMatrix N_m, B_m;
  FloatArray u_m;
 
    this->computeMicromorphicNMatrixAt(gp, N_m);
    this->computeMicromorphicBMatrixAt(gp, B_m);
    /// @todo generalization for general micromorphic continua -- should be parameter of this function
    this->giveElement()->computeVectorOf(IdMask_m, VM_Total, tStep, u_m);
    micromorphicVar.beProductOf(N_m, u_m);
    /*    micromorphicVar.at(1) += 1;
    micromorphicVar.at(2) += 1;
    micromorphicVar.at(3) += 1;
    */
    micromorphicVarGrad.beProductOf(B_m, u_m);
}


void 
QPlaneStrainGradLoc :: computeMicromorphicNMatrixAt(GaussPoint *gp,FloatMatrix &answer)
{
  this->computeMicromorphicBMatrixAt(gp, answer);
}


void
QPlaneStrainGradLoc :: computeMicromorphicBMatrixAt(GaussPoint *gp, FloatMatrix &answer)         
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix d2Ndx2; 
    interp->evald2Ndx2( d2Ndx2, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(8, this->giveNumberOfMicromorphicDofs());
    answer.zero();

    for ( int i = 1; i <= d2Ndx2.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = d2Ndx2.at(i, 1);
        answer.at(2, i * 2 - 1) = d2Ndx2.at(i, 3);

     	answer.at(3, i * 2 - 0) = d2Ndx2.at(i, 3);
        answer.at(4, i * 2 - 0) = d2Ndx2.at(i, 2);

        answer.at(5, i * 2 - 1) = d2Ndx2.at(i, 3);
        answer.at(6, i * 2 - 1) = d2Ndx2.at(i, 2);

        answer.at(7, i * 2 - 0) = d2Ndx2.at(i, 1);
        answer.at(8, i * 2 - 0) = d2Ndx2.at(i, 3);
    }
    
}
  
} // end namespace oofem
