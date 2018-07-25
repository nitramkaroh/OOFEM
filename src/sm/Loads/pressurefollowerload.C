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

namespace oofem {
REGISTER_BoundaryCondition(PressureFollowerLoad);

IRResultType PressureFollowerLoad :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    IR_GIVE_FIELD(ir, this->pressure, _IFT_PressureFollowerLoad_pressure);
    this->useTangent = ir->hasField(_IFT_PressureFollowerLoad_useTangent);
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

        e->giveInterpolation()->boundarySurfaceGiveNodes(bNodes, boundary);

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

        e->giveInterpolation()->boundarySurfaceGiveNodes(bNodes, boundary);

        e->giveBoundaryLocationArray(loc, bNodes, this->dofs, s, & masterdofids);
        this->computeLoadVectorFromElement(fe, e, boundary, tStep);
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
    IntegrationRule *iRule = pfli->surfaceGiveIntegrationRule(this->giveApproxOrder(), iSurf);

    IntArray bNodes;
    e->giveBoundarySurfaceNodes (bNodes, iSurf);
    double nNodes = bNodes.giveSize();

    for ( GaussPoint *gp : *iRule) {
        FloatMatrix dNdxi, dNk, dNe;    
        pfli->surfaceEvaldNdxi(dNdxi, iSurf, gp);
	this->giveSurfacedNdxi_dDdDMatrces(dNk, dNe, dNdxi, nNodes);
	// compute normal vector
        FloatArray n, dxdk,dxde;
	pfli->surfaceEvalDeformedNormalAt(n, dxdk, dxde, iSurf, gp, tStep);
	// compute surface N matirx
	FloatMatrix N;
	pfli->surfaceEvalNmatrixAt(N, iSurf, gp);
	FloatMatrix dXdk, dXde, dXedNk, dXkdNe;
	dXdk.giveMatrixOfAxialVector(dxdk);
	dXde.giveMatrixOfAxialVector(dxde);

	dXedNk.beProductOf(dXde,dNk);
	dXkdNe.beProductOf(dXdk,dNe);
	
	dXedNk.subtract(dXkdNe);
	

	double w = gp->giveWeight();	
	answer.plusProductUnsym(N, dXedNk, w);


    }
    answer.times(pressure);


}



void
PressureFollowerLoad :: giveSurfacedNdxi_dDdDMatrces(FloatMatrix &dNk_dD, FloatMatrix &dNe_dD, const FloatMatrix &dNdxi, int nNodes)
{

  //@todo nNodes*3 should be replaced by nDofs
  dNk_dD.resize(3,3*nNodes);
  dNe_dD.resize(3,3*nNodes);

  for( int i = 1; i <= nNodes; i++) {    
    dNk_dD.at(1,3*(i-1)+1) = dNdxi.at(i,1);
    dNk_dD.at(2,3*(i-1)+2) = dNdxi.at(i,1);
    dNk_dD.at(3,3*(i-1)+3) = dNdxi.at(i,1);

    dNe_dD.at(1,3*(i-1)+1) = dNdxi.at(i,2);
    dNe_dD.at(2,3*(i-1)+2) = dNdxi.at(i,2);
    dNe_dD.at(3,3*(i-1)+3) = dNdxi.at(i,2);
  }
    


}


  

void PressureFollowerLoad :: computeLoadVectorFromElement(FloatArray &answer, Element *e, int iSurf, TimeStep *tStep)
{

  PressureFollowerLoadElementInterface *pfli = static_cast< PressureFollowerLoadElementInterface * >(e->giveInterface(PressureFollowerLoadElementInterfaceType) );

    if ( !pfli ) {
        OOFEM_ERROR("Element doesn't implement the required Pressure Follower load interface!");
    }

    if ( iSurf == -1 ) {
      iSurf = 1;
    }

    answer.clear();
    IntegrationRule *iRule = pfli->surfaceGiveIntegrationRule(this->giveApproxOrder(), iSurf);
    double a = 0;
    for ( GaussPoint *gp : *iRule) {
	// compute normal vector
        FloatArray n, dxdksi,dxdeta;
	pfli->surfaceEvalDeformedNormalAt(n, dxdksi, dxdeta, iSurf, gp, tStep);
  	// compute surface N matirx
	FloatMatrix N;
	pfli->surfaceEvalNmatrixAt(N, iSurf, gp);
	// compute pressure follower load matrix
	double dV = pfli->surfaceEvalVolumeAround(gp, iSurf);
	//double w = gp->giveWeight();
	answer.plusProduct(N, n, dV);
	a += dV*n.computeNorm();
	
	
    }
    // ask time distribution
    double factor = this->giveTimeFunction()->evaluate(tStep, VM_Total);
    answer.times(pressure*factor);

    
}


 

} // end namespace oofem
