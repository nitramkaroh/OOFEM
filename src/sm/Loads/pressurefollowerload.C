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

#include "pressurefollowerload.h"
#include "element.h"
#include "node.h"
#include "crosssection.h"
#include "error.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "classfactory.h"
#include "set.h"
#include "sparsemtrx.h"
#include "function.h"


#include "Elements/structural3delement.h"



#include "integrationdomain.h"
#include "integrationrule.h"
#include "gausspoint.h"
#include "mathfem.h"

#include <utility>
#include <list>
#include <memory>


#include "../sm/Elements/Axisymmetry/l4axisymm.h"

namespace oofem {
REGISTER_BoundaryCondition(PressureFollowerLoad);

IRResultType PressureFollowerLoad :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    IR_GIVE_FIELD(ir, this->pressure, _IFT_PressureFollowerLoad_pressure);
    this->nfl = ir->hasField(_IFT_PressureFollowerLoad_nfl);
    if(nfl == false) {
      this->useTangent = ir->hasField(_IFT_PressureFollowerLoad_useTangent);
    } else {
      this->useTangent = false;
    }
    
    return ActiveBoundaryCondition :: initializeFrom(ir);
}

void PressureFollowerLoad :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                           const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
        return;
    }

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    IntArray bNodes;

    rows.resize(boundaries.giveSize() / 2);
    cols.resize(boundaries.giveSize() / 2);

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        e->giveBoundaryLocationArray(rows [ pos ], bNodes, this->dofs, r_s);
        e->giveBoundaryLocationArray(cols [ pos ], bNodes, this->dofs, c_s);
    }
}

void PressureFollowerLoad :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                                 CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( !this->useTangent || type != TangentStiffnessMatrix ) {
        return;
    }


    FloatMatrix Ke;
    IntArray r_loc, c_loc, bNodes;

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
	if(e->giveSpatialDimension() == 2) {
	  e->giveInterpolation()->boundaryEdgeGiveNodes(bNodes, boundary);
	} else {
	  e->giveInterpolation()->boundarySurfaceGiveNodes(bNodes, boundary);
	}

        e->giveBoundaryLocationArray(r_loc, bNodes, this->dofs, r_s);
        e->giveBoundaryLocationArray(c_loc, bNodes, this->dofs, c_s);
        this->computeTangentFromElement(Ke, e, boundary, tStep);
        answer.assemble(r_loc, c_loc, Ke);
    }
}

void PressureFollowerLoad :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                       CharType type, ValueModeType mode,
                                                       const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    if ( type != ExternalForcesVector ) {
        return;
    }

    FloatArray fe;
    IntArray loc, masterdofids, bNodes;

    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
	if(e->giveSpatialDimension() == 2) {
	  e->giveInterpolation()->boundaryEdgeGiveNodes(bNodes, boundary);
	} else {
	  e->giveInterpolation()->boundarySurfaceGiveNodes(bNodes, boundary);
	}
        e->giveBoundaryLocationArray(loc, bNodes, this->dofs, s, & masterdofids);
        this->computeLoadVectorFromElement(fe, e, boundary, tStep, mode);
        answer.assemble(fe, loc);
        if ( eNorms ) {
            eNorms->assembleSquared(fe, masterdofids);
        }
    }
}

void PressureFollowerLoad :: computeTangentFromElement(FloatMatrix &answer, Element *e, int iSurf, TimeStep *tStep)
{
  // compute pressure follower load matrix
    PressureFollowerLoadElementInterface *pfli = static_cast< PressureFollowerLoadElementInterface * >(e->giveInterface(PressureFollowerLoadElementInterfaceType) );

    if ( !pfli ) {
        OOFEM_ERROR("Element doesn't implement the required Pressure Follower load interface!");
    }

    if ( iSurf == -1 ) {
      iSurf = 1;
    }

    answer.clear();
    int order = 4;
    IntegrationRule *iRule = pfli->surfaceGiveIntegrationRule(order, iSurf);

    IntArray bNodes;
    e->giveBoundarySurfaceNodes(bNodes, iSurf);
    double nNodes = bNodes.giveSize();
    FloatMatrix K;
    FloatMatrix testAnswer(12,12);
    if( e->giveSpatialDimension() == 3) {
      e->giveBoundarySurfaceNodes (bNodes, iSurf);
      for ( GaussPoint *gp : *iRule) {
        FloatMatrix dN, dNksi, dNeta;    
        pfli->surfaceEvaldNdxi(dN, iSurf, gp);
	this->giveSurface_dNdKsi_dNdEta(dNksi, dNeta, dN, nNodes);
	// compute normal vector
        FloatArray n, dxdk,dxde;
	pfli->surfaceEvalDeformedNormalAt(n, dxdk, dxde, iSurf, gp, tStep);
	// compute surface N matirx
	FloatMatrix N;
	pfli->surfaceEvalNmatrixAt(N, iSurf, gp);
	FloatMatrix dXdk, dXde, dXedNk, dXkdNe;
	dXdk.giveMatrixOfAxialVector(dxdk);
	dXde.giveMatrixOfAxialVector(dxde);
	
	dXedNk.beProductOf(dXde,dNksi);
	dXkdNe.beProductOf(dXdk,dNeta);
	
	dXkdNe.subtract(dXedNk);
	
	
	double w = gp->giveWeight();	
	answer.plusProductUnsym(N, dXkdNe, w);

	pfli->surfaceEvalNumericalStiffMatrixAt(K, dxde, dxdk, 1, gp, tStep);
	testAnswer.add(K);
	//answer.add(K);
      }
    } else if ( e->giveSpatialDimension() == 2) {
      //  if(edge) {
	e->giveBoundaryEdgeNodes (bNodes, iSurf);
	for ( GaussPoint *gp : *iRule) { 
	  // compute surface/edge N matirx
	  FloatMatrix N, dN, dNI;
	  FloatMatrix eps2 = {{0, -1},{1,0}};
	  pfli->surfaceEvalNmatrixAt(N, iSurf, gp);
	  pfli->surfaceEvaldNdxi(dN, iSurf, gp);
	  FloatMatrix eNksi;
	  eNksi.beProductOf(eps2,dN);
	  double w = gp->giveWeight();
	  answer.plusProductUnsym(N, eNksi, w);
	  FloatArray n, dxdk,dxde;
	  FloatMatrix K, test, dNdx;
	  test.beTProductOf(N, eNksi);
	  pfli->surfaceEvalDeformedNormalAt(n, dxdk, dxde, iSurf, gp, tStep);
	  pfli->surfaceEvalNumericalStiffMatrixAt(K, dNdx, dxde, dxdk, 1, gp, tStep);

	}
	
	// }



	/*} else { // surface in 2d, i.e., membrane???
	e->giveBoundaryEdgeNodes (bNodes, iSurf);
	for ( GaussPoint *gp : *iRule) { 
	  // compute surface/edge N matirx
	  FloatMatrix N, dndu;
	  pfli->surfaceEvalNmatrixAt(N, iSurf, gp);
	  pfli->surfaceEvalNormalDerivative(dndu, iSurf, gp, tStep);
	  double w = gp->giveWeight();	
	  answer.plusProductUnsym(N, dndu, w);
	}
	}*/
    }
    //answer = testAnswer;
    double factor = this->giveTimeFunction()->evaluate(tStep, VM_Total);
    // minus in 3d
    //answer.times(-pressure*factor);
    // plus in 2d ?
    answer.times(pressure*factor);
    /*L4Axisymm *le = dynamic_cast<L4Axisymm* > (e);
    e->giveInterpolation()->boundaryEdgeGiveNodes(bNodes, iSurf);
    FloatArray(vU);
    FloatMatrix test(4,4);
    le->computeBoundaryVectorOf(bNodes, {D_u, D_v}, VM_Total, tStep, vU); // solution vector    
    double b11 = 2./3.;
    double b12 = 1./3.;
    double b21 = 1./3.;
    double b22 = 2./3.;
    double r1 =  le->giveCellGeometryWrapper(tStep)->giveVertexCoordinates(bNodes.at(1))->at(1) + vU.at(1);
    double r2 =  le->giveCellGeometryWrapper(tStep)->giveVertexCoordinates(bNodes.at(2))->at(1) + vU.at(3);
    double z1 =  le->giveCellGeometryWrapper(tStep)->giveVertexCoordinates(bNodes.at(1))->at(2) + vU.at(2);
    double z2 =  le->giveCellGeometryWrapper(tStep)->giveVertexCoordinates(bNodes.at(2))->at(2) + vU.at(4);

    double a1 = b11*r1+b12*r2;
    double a2 = b21*r1+b22*r2;

    
    test = {{b11*(z2-z1), -b11*(r2-r1) - a1, b12*(z2-z1), -b12*(r2-r1)-a2},{a1, 0, a2, 0},{b12*(z2-z1),-b12*(r2-r1)+a1, b22*(z2-z1),-b22*(r2-r1)+a2},{-a1, 0, -a2, 0}};
    test.times(-0.5*pressure);

    */

}



void
PressureFollowerLoad :: giveSurface_dNdKsi_dNdEta(FloatMatrix &dNdksi, FloatMatrix &dNdeta, const FloatMatrix &dN, int nNodes)
{

  //@todo nNodes*3 should be replaced by nDofs
  dNdksi.resize(3,3*nNodes);
  dNdeta.resize(3,3*nNodes);

  for( int i = 1; i <= nNodes; i++) {    
    dNdksi.at(1,3*(i-1)+1) = dN.at(i,1);
    dNdksi.at(2,3*(i-1)+2) = dN.at(i,1);
    dNdksi.at(3,3*(i-1)+3) = dN.at(i,1);

    dNdeta.at(1,3*(i-1)+1) = dN.at(i,2);
    dNdeta.at(2,3*(i-1)+2) = dN.at(i,2);
    dNdeta.at(3,3*(i-1)+3) = dN.at(i,2);
  }
    


}


  

  void PressureFollowerLoad :: computeLoadVectorFromElement(FloatArray &answer, Element *e, int iSurf, TimeStep *tStep, ValueModeType mode)
{

  PressureFollowerLoadElementInterface *pfli = static_cast< PressureFollowerLoadElementInterface * >(e->giveInterface(PressureFollowerLoadElementInterfaceType) );

    if ( !pfli ) {
        OOFEM_ERROR("Element doesn't implement the required Pressure Follower load interface!");
    }

    if ( iSurf == -1 ) {
      iSurf = 1;
    }

    answer.clear();
    int order = 4;
    IntegrationRule *iRule = pfli->surfaceGiveIntegrationRule(order, iSurf);
    double a = 0;
    for ( GaussPoint *gp : *iRule) {
	// compute normal vector
        FloatArray n, dxdksi,dxdeta;
	if(nfl == false) {
	  pfli->surfaceEvalDeformedNormalAt(n, dxdksi, dxdeta, iSurf, gp, tStep);
	} else {
	  pfli->surfaceEvalNormalAt(n, dxdksi, dxdeta, iSurf, gp, tStep);
	}
  	// compute surface N matirx
	FloatMatrix N;
	pfli->surfaceEvalNmatrixAt(N, iSurf, gp);
	// compute pressure follower load matrix
	//double dV = pfli->surfaceEvalVolumeAround(gp, iSurf);
	// @todo parametrization [-1 1] => factor 2 ???
	double w = gp->giveWeight();
	answer.plusProduct(N, n, w);
	a += n.computeNorm();
	
	
    }
    // ask time distribution
    double factor = this->giveTimeFunction()->evaluate(tStep, VM_Total);
    answer.times(-pressure*factor);


    /*
    L4Axisymm *le = dynamic_cast<L4Axisymm* > (e);

    FloatArray(vU);
    FloatArray test(4);

    IntArray bNodes;
    
    e->giveInterpolation()->boundaryEdgeGiveNodes(bNodes, iSurf);


    
    le->computeBoundaryVectorOf(bNodes, {D_u, D_v}, VM_Total, tStep, vU); // solution vector    
    double b11 = 2./3.;
    double b12 = 1./3.;
    double b21 = 1./3.;
    double b22 = 2./3.;
    double r1 =  le->giveCellGeometryWrapper(tStep)->giveVertexCoordinates(bNodes.at(1))->at(1) + vU.at(1);
    double r2 =  le->giveCellGeometryWrapper(tStep)->giveVertexCoordinates(bNodes.at(2))->at(1) + vU.at(3);
    double z1 =  le->giveCellGeometryWrapper(tStep)->giveVertexCoordinates(bNodes.at(1))->at(2) + vU.at(2);
    double z2 =  le->giveCellGeometryWrapper(tStep)->giveVertexCoordinates(bNodes.at(2))->at(2) + vU.at(4);

    double a1 = b11*r1+b12*r2;
    double a2 = b21*r1+b22*r2;

    
    test = {a1*(z2-z1), -a1*(r2-r1), a2*(z2-z1), -a2*(r2-r1)};
    test.times(-0.5*pressure);
    */    
}


 

} // end namespace oofem
