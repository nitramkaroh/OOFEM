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

#include "../sm/Elements/PlaneStress/quadmembrane.h"
#include "../sm/Materials/structuralms.h"
#include "fei2dquadlin.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "strainvector.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(QuadMembrane);

FEI2dQuadLin QuadMembrane :: interpolation(1, 2);

QuadMembrane :: QuadMembrane(int n, Domain *aDomain) :
    PlaneStress2d(n, aDomain)
    //, PressureFollowerLoadElementInterface(this)
    // Constructor.
{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 4;
    nlGeometry = 1;
}

QuadMembrane :: ~QuadMembrane()
// Destructor
{ }

FEInterpolation *QuadMembrane :: giveInterpolation() const { return & interpolation; }

void
QuadMembrane :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, int li, int ui)
//
// Returns the [3x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );

    answer.resize(3, 8);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);
    }

  
  OOFEM_ERROR("Should be called???")
}


void
QuadMembrane :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, double alpha)
//
// Returns the [4x8] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// @todo not checked if correct
{
    FloatMatrix dNdx;

    this->interpolation.evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
    answer.resize(6, 12);

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);     // du/dx - 1
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);     // dv/dy - 2
	answer.at(3, 3 * i - 2) = dNdx.at(i, 2);     // du/dy - 6
	answer.at(4, 3 * i - 0) = dNdx.at(i, 2);     // dw/dy - 7
        answer.at(5, 3 * i - 0) = dNdx.at(i, 1);     // dw/dx - 8
	answer.at(6, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx - 9
    }

}



void
QuadMembrane :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
      D_u, D_v, D_w
    };
}

  

IRResultType
QuadMembrane :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 4;
    IRResultType result = PlaneStressElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    if ( numberOfGaussPoints != 1 && numberOfGaussPoints != 4 && numberOfGaussPoints != 9 && numberOfGaussPoints != 16 && numberOfGaussPoints != 25 ) {
        numberOfGaussPoints = 4;
        OOFEM_WARNING("Number of Gauss points enforced to 4");
    }

    return IRRT_OK;
}


void
QuadMembrane :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType modeType)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    this->computeVectorOf({D_u, D_v, D_w}, modeType, tStep, u); // solution vector
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // Displacement gradient H = du/dX
    FloatArray rF, n, a, b;
    FloatMatrix B;
    
    this->computeBHmatrixAt(gp, B, tStep, 0);
    rF.beProductOf(B, u);
    /* this->surfaceEvalDeformedNormalAt(n, a, b, 1, gp, tStep);
    if(n.at(3) < 0) {
      n.times(-1.0);
    }
    
    n.normalize();
    answer.resize(9);
    
    answer.at(1) = rF.at(1) + 1.0;
    answer.at(2) = rF.at(2) + 1.0;
    answer.at(3) = n.at(3);
    answer.at(4) = n.at(2);
    answer.at(5) = n.at(1);
    answer.at(6) = rF.at(3);
    answer.at(7) = rF.at(4);
    answer.at(8) = rF.at(5);
    answer.at(9) = rF.at(6);
    */
    answer.resize(6);
    answer.at(1) = rF.at(1) + 1.0;
    answer.at(2) = rF.at(2) + 1.0;
    answer.at(3) = rF.at(3);
    answer.at(4) = rF.at(4);
    answer.at(5) = rF.at(5);
    answer.at(6) = rF.at(6);

    
}






// support for pressure follower load interface
/*Interface*
QuadMembrane :: giveInterface(InterfaceType interface)
{
    if ( interface == PressureFollowerLoadElementInterfaceType) {
        return static_cast< PressureFollowerLoadElementInterface* >(this);
    }
    return NULL;
}
*/

void
QuadMembrane :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    int numNodes = this->giveNumberOfDofManagers();
    FloatArray N(numNodes);

    int dim = 3;

    answer.resize(dim, dim * numNodes);
    answer.zero();
    giveInterpolation()->evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );

    answer.beNMatrixOf(N, dim);
}

  
void
QuadMembrane ::  surfaceEvaldNdxi(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{

  FEInterpolation2d *interp2d = static_cast<FEInterpolation2d*> (this->giveInterpolation());
  interp2d->boundarySurfaceEvaldNdx(answer, iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));
}


void
QuadMembrane ::  surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdksi, int iSurf, GaussPoint *gp, TimeStep *tStep)
{
  IntArray bNodes;
  FloatArray gcoords, lcoords, vU;
  FloatMatrix dNdxi, dxdxi, x;

  lcoords = gp->giveNaturalCoordinates();
  this->giveBoundarySurfaceNodes (bNodes, iSurf);

  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);
  double nNodes = bNodes.giveSize();
  x.resize(nNodes,3);
  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    // compute actual node positions for Total Lagrangean formulation
    this->computeBoundaryVectorOf(bNodes, {D_u, D_v, D_w}, VM_Total, tStep, vU); // solution vector    
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(i,1) = node->giveCoordinate(1) + vU.at( (i-1) * 3 + 1);
      x.at(i,2) = node->giveCoordinate(2) + vU.at( (i-1) * 3 + 2);
      x.at(i,3) = node->giveCoordinate(3) + vU.at( (i-1) * 3 + 3);
    }
  } else {
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(i,1) = node->giveCoordinate(1);
      x.at(i,2) = node->giveCoordinate(2);
      x.at(i,3) = node->giveCoordinate(3);
    }
    
  }

  /* testing stuff
  FloatMatrix dNdk(4,1), xx, dxdx;
  for(int i = 1; i <= 4; i++) {
    dNdk.at(i,1) = dNdxi.at(i,1);
  }

  dxdx.beTProductOf(x, dNdk);
  */
  
  dxdxi.beTProductOf(dNdxi,x);
  dxdxi.copyRow(dxdksi, 2);
  dxdxi.copyRow(dxdeta, 1);
  answer.beVectorProductOf(dxdksi, dxdeta);

}
  


} // end namespace oofem
