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

#include "../sm/Elements/Micromorphic/GradientPolyconvex/planestraingradpolyconvex.h"
#include "fei2dquadlin.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(PlaneStrainGradPolyconvex);

FEI2dQuadLin PlaneStrainGradPolyconvex :: interpolation(1, 2);

PlaneStrainGradPolyconvex :: PlaneStrainGradPolyconvex(int n, Domain *aDomain) : Quad1PlaneStrain(n, aDomain), BaseMicromorphicElement()
    // Constructor.
{

  displacementDofsOrdering = {1,2,8,9,15,16,22,23};
  micromorphicDofsOrdering = {3,4,5,6,7,10,11,12,13,14,17,18,19,20,21,24,25,26,27,28};
}


void 
PlaneStrainGradPolyconvex :: computeMicromorphicNMatrixAt(GaussPoint *gp,FloatMatrix &answer)
{
    FloatArray n;
    this->interpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 5);
}


void
PlaneStrainGradPolyconvex :: computeMicromorphicBMatrixAt(GaussPoint *gp, FloatMatrix &answer)         
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(10, this->giveNumberOfMicromorphicDofs());
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {

        answer.at(1, i * 5 - 4) = dNdx.at(i, 1);
        answer.at(2, i * 5 - 4) = dNdx.at(i, 2);

     	answer.at(3, i * 5 - 3) = dNdx.at(i, 1);
        answer.at(4, i * 5 - 3) = dNdx.at(i, 2);

        answer.at(5, i * 5 - 2) = dNdx.at(i, 1);
        answer.at(6, i * 5 - 2) = dNdx.at(i, 2);

        answer.at(7, i * 5 - 1) = dNdx.at(i, 1);
        answer.at(8, i * 5 - 1) = dNdx.at(i, 2);
	
        answer.at(9, i * 5 - 0) = dNdx.at(i, 1);
        answer.at(10, i * 5 - 0) = dNdx.at(i, 2);
    }
    
}


void
PlaneStrainGradPolyconvex :: computeMicromorphicVars(FloatArray &micromorphicVar,FloatArray &micromorphicVarGrad, IntArray IdMask_m, GaussPoint *gp, TimeStep *tStep)
{
  FloatMatrix N_m, B_m;
  FloatArray u_m;
 
    this->computeMicromorphicNMatrixAt(gp, N_m);
    this->computeMicromorphicBMatrixAt(gp, B_m);
    /// @todo generalization for general micromorphic continua -- should be parameter of this function
    this->giveElement()->computeVectorOf(IdMask_m, VM_Total, tStep, u_m);
    micromorphicVar.beProductOf(N_m, u_m);
    micromorphicVar.at(1) += 1;
    micromorphicVar.at(2) += 1;
    micromorphicVar.at(3) += 1;
    micromorphicVarGrad.beProductOf(B_m, u_m);
}
 


void
PlaneStrainGradPolyconvex :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {D_u, D_v, M_X11, M_X22, M_X33, M_X12, M_X21};
}


void
PlaneStrainGradPolyconvex :: giveDofManDofIDMask_u(IntArray &answer)
{
  answer = {D_u, D_v};
}


void
PlaneStrainGradPolyconvex :: giveDofManDofIDMask_m(IntArray &answer)
{
  answer = {M_X11, M_X22, M_X33, M_X12, M_X21};
}


void
PlaneStrainGradPolyconvex ::  postInitialize() 
{
  BaseMicromorphicElement :: postInitialize();
  Quad1PlaneStrain :: postInitialize();
}

}
