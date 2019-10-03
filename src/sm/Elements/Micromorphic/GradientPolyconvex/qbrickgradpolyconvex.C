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

#include "../sm/Elements/Micromorphic/GradientPolyconvex/qbrickgradpolyconvex.h"
#include "fei3dhexaquad.h"

#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(QBrickGradPolyconvex);
  
FEI3dHexaQuad QBrickGradPolyconvex :: q_interpolation;
FEI3dHexaLin  QBrickGradPolyconvex :: l_interpolation;

QBrickGradPolyconvex :: QBrickGradPolyconvex(int n, Domain *aDomain) : QSpace(n, aDomain), BaseMicromorphicElement()
    // Constructor.
{
  int index = 0;
  for(int iNode = 1; iNode <= 8; iNode++) {
    for( int iDof = 1; iDof <= 12; iDof++ ) {
      index++;
      if(iDof <= 3) {
	displacementDofsOrdering.followedBy(index);
      } else {
	micromorphicDofsOrdering.followedBy(index);
      }
	  
    }

  }


  
  for(int iNode = 9; iNode <= 20; iNode++) {
    for( int iDof = 1; iDof <= 3; iDof++ ) {
      index++;
      displacementDofsOrdering.followedBy(index);
    }
  }


  
  /*  displacementDofsOrdering = {1,2,3,13,14,15,8,9,15,16,22,23};
  micromorphicDofsOrdering = {3,4,5,6,7,10,11,12,13,14,17,18,19,20,21,24,25,26,27,28};
  */
}


void 
LBrickGradPolyconvex :: computeMicromorphicNMatrixAt(GaussPoint *gp,FloatMatrix &answer)
{
    FloatArray n;
    this->interpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 9);
}


void
LBrickGradPolyconvex :: computeMicromorphicBMatrixAt(GaussPoint *gp, FloatMatrix &answer)         
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(27, this->giveNumberOfMicromorphicDofs());
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {

        answer.at(1, i * 9 - 8) = dNdx.at(i, 1);
        answer.at(2, i * 9 - 8) = dNdx.at(i, 2);
	answer.at(3, i * 9 - 8) = dNdx.at(i, 3);


	answer.at(4, i * 9 - 7) = dNdx.at(i, 1);
        answer.at(5, i * 9 - 7) = dNdx.at(i, 2);
	answer.at(6, i * 9 - 7) = dNdx.at(i, 3);

	answer.at(7, i * 9 - 6) = dNdx.at(i, 1);
        answer.at(8, i * 9 - 6) = dNdx.at(i, 2);
	answer.at(9, i * 9 - 6) = dNdx.at(i, 3);

	answer.at(10, i * 9 - 5) = dNdx.at(i, 1);
        answer.at(11, i * 9 - 5) = dNdx.at(i, 2);
	answer.at(12, i * 9 - 5) = dNdx.at(i, 3);

	answer.at(13, i * 9 - 4) = dNdx.at(i, 1);
        answer.at(14, i * 9 - 4) = dNdx.at(i, 2);
	answer.at(15, i * 9 - 4) = dNdx.at(i, 3);

	answer.at(16, i * 9 - 3) = dNdx.at(i, 1);
        answer.at(17, i * 9 - 3) = dNdx.at(i, 2);
	answer.at(18, i * 9 - 3) = dNdx.at(i, 3);

	answer.at(19, i * 9 - 2) = dNdx.at(i, 1);
        answer.at(20, i * 9 - 2) = dNdx.at(i, 2);
	answer.at(21, i * 9 - 2) = dNdx.at(i, 3);

	answer.at(22, i * 9 - 1) = dNdx.at(i, 1);
        answer.at(23, i * 9 - 1) = dNdx.at(i, 2);
	answer.at(24, i * 9 - 1) = dNdx.at(i, 3);


	answer.at(25, i * 9 - 0) = dNdx.at(i, 1);
        answer.at(26, i * 9 - 0) = dNdx.at(i, 2);
	answer.at(27, i * 9 - 0) = dNdx.at(i, 3);
	
    }
    
}


void
LBrickGradPolyconvex :: computeMicromorphicVars(FloatArray &micromorphicVar,FloatArray &micromorphicVarGrad, IntArray IdMask_m, GaussPoint *gp, TimeStep *tStep)
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
LBrickGradPolyconvex :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {D_u, D_v, D_w, M_X11, M_X22, M_X33, M_X23, M_X13, M_X12, M_X32, M_X31, M_X21};
}


void
LBrickGradPolyconvex :: giveDofManDofIDMask_u(IntArray &answer)
{
  answer = {D_u, D_v, D_w};
}


void
LBrickGradPolyconvex :: giveDofManDofIDMask_m(IntArray &answer)
{

  answer = {M_X11, M_X22, M_X33, M_X23, M_X13, M_X12, M_X32, M_X31, M_X21};

}


void
LBrickGradPolyconvex ::  postInitialize() 
{
  BaseMicromorphicElement :: postInitialize();
  LSpace :: postInitialize();
}

}












#include "Elements/3D/qspace.h"
#include "CrossSections/structuralcrosssection.h"
#include "fei3dhexaquad.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(QSpace);

FEI3dHexaQuad QSpace :: interpolation;

QSpace :: QSpace(int n, Domain *aDomain) : Structural3DElement(n, aDomain), ZZNodalRecoveryModelInterface(this)
{
    numberOfDofMans = 20;
}


IRResultType
QSpace :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 27;

    return Structural3DElement :: initializeFrom(ir);
}


FEInterpolation *QSpace :: giveInterpolation() const { return & interpolation; }

// ******************************
// ***  Surface load support  ***
// ******************************

IntegrationRule *
QSpace :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->SetUpPointsOnSquare(npoints, _Unknown);
    return iRule;
}


int
QSpace :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{
    // returns transformation matrix from
    // surface local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)

    // definition of local c.s on surface:
    // local z axis - perpendicular to surface, pointing outwards from element
    // local x axis - is in global xy plane (perpendicular to global z axis)
    // local y axis - completes the righ hand side cs.

    /*
     * OOFEM_ERROR("surface local coordinate system not supported");
     * return 1;
     */
    FloatArray gc(3);
    FloatArray h1(3), h2(3), nn(3), n(3);
    IntArray snodes(4);

    answer.resize(3, 3);
    answer.zero();

    this->interpolation.computeSurfaceMapping(snodes, dofManArray, iSurf);
    for ( int i = 1; i <= 4; i++ ) {
        gc.add( * domain->giveNode( snodes.at(i) )->giveCoordinates() );
    }

    gc.times(1. / 4.);
    // determine "average normal"
    for ( int i = 1; i <= 4; i++ ) {
        int j = ( i ) % 4 + 1;
        h1.beDifferenceOf(* domain->giveNode( snodes.at(i) )->giveCoordinates(), gc);
        h2.beDifferenceOf(* domain->giveNode( snodes.at(j) )->giveCoordinates(), gc);
        n.beVectorProductOf(h1, h2);
        if ( n.computeSquaredNorm() > 1.e-6 ) {
            n.normalize();
        }

        nn.add(n);
    }

    nn.times(1. / 4.);
    if ( nn.computeSquaredNorm() < 1.e-6 ) {
        answer.zero();
    }

    nn.normalize();
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, 3) = nn.at(i);
    }

    // determine lcs of surface
    // local x axis in xy plane
    double test = fabs(fabs( nn.at(3) ) - 1.0);
    if ( test < 1.e-5 ) {
        h1.at(1) = answer.at(1, 1) = 1.0;
        h1.at(2) = answer.at(2, 1) = 0.0;
    } else {
        h1.at(1) = answer.at(1, 1) = answer.at(2, 3);
        h1.at(2) = answer.at(2, 1) = -answer.at(1, 3);
    }

    h1.at(3) = answer.at(3, 1) = 0.0;
    // local y axis perpendicular to local x,z axes
    h2.beVectorProductOf(nn, h1);
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, 2) = h2.at(i);
    }

    return 1;
}

Interface *
QSpace :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    OOFEM_LOG_INFO("Interface on Qspace element not supported");
    return NULL;
}

void
QSpace :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(20);
    for ( int i = 1; i <= 20; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
QSpace :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 20; i++ ) {
        if ( this->giveNode(i)->giveNumber() == pap ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("unknown node number %d", pap);
    }
}

int
QSpace :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints();
}


SPRPatchType
QSpace :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiQuadratic;
}


void
QSpace :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    answer.clear();
    OOFEM_WARNING("IP values will not be transferred to nodes. Use ZZNodalRecovery instead (parameter stype 1)");
}

} // end namespace oofem
