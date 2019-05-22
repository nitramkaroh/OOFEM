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

#include "../sm/Elements/Micromorphic/Micropolar/lbrickmicropolar.h"
#include "fei3dhexalin.h"

#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(LBrickMicropolar);
  
FEI3dHexaLin LBrickMicropolar :: interpolation;

LBrickMicropolar :: LBrickMicropolar(int n, Domain *aDomain) : LSpace(n, aDomain), BaseMicromorphicElement()
    // Constructor.
{
  int index = 0;
  for(int iNode = 1; iNode <=8; iNode++) {
    for( int iDof = 1; iDof <= 6; iDof++ ) {
      index++;
      if(iDof <= 3) {
	displacementDofsOrdering.followedBy(index);
      } else {
	micromorphicDofsOrdering.followedBy(index);
      }
	  
    }

  }
  /*  displacementDofsOrdering = {1,2,3,13,14,15,8,9,15,16,22,23};
  micromorphicDofsOrdering = {3,4,5,6,7,10,11,12,13,14,17,18,19,20,21,24,25,26,27,28};
  */
}


void 
LBrickMicropolar :: computeMicromorphicNMatrixAt(GaussPoint *gp,FloatMatrix &answer)
{
    FloatArray n;
    this->interpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 3);
}


void
LBrickMicropolar :: computeMicromorphicBMatrixAt(GaussPoint *gp, FloatMatrix &answer)         
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(9, this->giveNumberOfMicromorphicDofs());
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {

	answer.at(1, i * 3 - 2) = dNdx.at(i, 1);
        answer.at(2, i * 3 - 2) = dNdx.at(i, 2);
	answer.at(3, i * 3 - 2) = dNdx.at(i, 3);

	answer.at(4, i * 3 - 1) = dNdx.at(i, 1);
        answer.at(5, i * 3 - 1) = dNdx.at(i, 2);
	answer.at(6, i * 3 - 1) = dNdx.at(i, 3);


	answer.at(7, i * 3 - 0) = dNdx.at(i, 1);
        answer.at(8, i * 3 - 0) = dNdx.at(i, 2);
	answer.at(9, i * 3 - 0) = dNdx.at(i, 3);
	
    }
    
}


void
LBrickMicropolar :: computeMicromorphicVars(FloatArray &micromorphicVar,FloatArray &micromorphicVarGrad, IntArray IdMask_m, GaussPoint *gp, TimeStep *tStep)
{
  FloatMatrix N_m, B_m;
  FloatArray u_m;
 
    this->computeMicromorphicNMatrixAt(gp, N_m);
    this->computeMicromorphicBMatrixAt(gp, B_m);
    /// @todo generalization for general micromorphic continua -- should be parameter of this function
    this->giveElement()->computeVectorOf(IdMask_m, VM_Total, tStep, u_m);
    micromorphicVar.beProductOf(N_m, u_m);
    micromorphicVarGrad.beProductOf(B_m, u_m);
}
 


void
LBrickMicropolar :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {D_u, D_v, D_w, M_W1, M_W2, M_W3};
}


void
LBrickMicropolar :: giveDofManDofIDMask_u(IntArray &answer)
{
  answer = {D_u, D_v, D_w};
}


void
LBrickMicropolar :: giveDofManDofIDMask_m(IntArray &answer)
{

  answer = {M_W1, M_W2, M_W3};

}


void
LBrickMicropolar ::  postInitialize() 
{
  BaseMicromorphicElement :: postInitialize();
  LSpace :: postInitialize();
}

}


