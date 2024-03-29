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

#include "Elements/structural2delement.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "CrossSections/structuralcrosssection.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"

namespace oofem {
Structural2DElement :: Structural2DElement(int n, Domain *aDomain) :
  NLStructuralElement(n, aDomain), FbarElementExtensionInterface(aDomain),PressureFollowerLoadElementInterface(this), matRotation(false)
{
    cellGeometryWrapper = NULL;
}

Structural2DElement :: ~Structural2DElement()
{
    if ( cellGeometryWrapper ) {
        delete cellGeometryWrapper;
    }
}


IRResultType
Structural2DElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result = NLStructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    
    result = FbarElementExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
      return result;
    }

    matRotation = ir->hasField(_IFT_Structural2DElement_materialCoordinateSystem); //|| this->elemLocalCS.isNotEmpty();
    return IRRT_OK;
}


Interface *
Structural2DElement :: giveInterface(InterfaceType interface)
{
   if ( interface == PressureFollowerLoadElementInterfaceType) {
        return static_cast< PressureFollowerLoadElementInterface* >(this);
    }   
    return NLStructuralElement::giveInterface(interface);
}

  
void
Structural2DElement :: postInitialize()
{
    // Element must be created before giveNumberOfNodes can be called
    StructuralElement :: postInitialize();
    this->numberOfDofMans = this->giveNumberOfNodes();
}


int
Structural2DElement :: giveNumberOfNodes() const
{
    return this->giveInterpolation()->giveNumberOfNodes();
}


FEICellGeometry *
Structural2DElement :: giveCellGeometryWrapper(TimeStep *tStep, double alpha)
{
    if ( !cellGeometryWrapper ) {
      if(this->giveDomain()->giveEngngModel()->giveFormulation() == AL) {
	cellGeometryWrapper = new FEIElementDeformedGeometryWrapper(this);
      } else {
        cellGeometryWrapper = new FEIElementGeometryWrapper(this);
      }
    }

    if(this->giveDomain()->giveEngngModel()->giveFormulation() == AL) {
      FEIElementDeformedGeometryWrapper *dgw = static_cast<FEIElementDeformedGeometryWrapper*> (cellGeometryWrapper);      
      if(tStep != NULL) {
	dgw->setTimeStep(tStep);
      }
      dgw->setAlpha(alpha);
    }

    return cellGeometryWrapper;
}


int
Structural2DElement :: computeNumberOfDofs()
{
    ///@todo move one hiearchy up and generalize
    IntArray dofIdMask;
    this->giveDofManDofIDMask(-1, dofIdMask); // ok for standard elements
    return this->giveInterpolation()->giveNumberOfNodes() * dofIdMask.giveSize();
}


void
Structural2DElement :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v
    };
}


double
Structural2DElement :: computeVolumeAround(GaussPoint *gp)
{
    // Computes the volume element dV associated with the given gp.

    double weight = gp->giveWeight();
    const FloatArray &lCoords = gp->giveNaturalCoordinates(); // local/natural coords of the gp (parent domain)
    double detJ = fabs( this->giveInterpolation()->giveTransformationJacobian( lCoords, * this->giveCellGeometryWrapper() ) );
    double thickness = this->giveCrossSection()->give(CS_Thickness, gp); // the cross section keeps track of the thickness

    return detJ * thickness * weight; // dV
}


void
Structural2DElement :: computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    // Default: create one integration rule
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}


void
Structural2DElement :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType modeType)
{
  if(!FbarFlag) {
    NLStructuralElement :: computeDeformationGradientVector(answer, gp, tStep, modeType);
  } else {
    FbarElementExtensionInterface :: computeFbarDeformationGradientVector(answer, gp, tStep, this, modeType);
  }
}

void
Structural2DElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  if(!FbarFlag) {
    NLStructuralElement :: computeStiffnessMatrix(answer, rMode, tStep);
  } else {
    FbarElementExtensionInterface :: computeFbarStiffnessMatrix(answer, rMode, tStep, this);
  }
}

void
Structural2DElement :: computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the first Piola-Kirchhoff stress vector containing the stresses at the  Gauss point
    // gp of the receiver at time step tStep. The deformation gradient F is computed and sent as
    // input to the material model.
    StructuralCrossSection *cs = this->giveStructuralCrossSection();

    FloatArray vF;
      if(!FbarFlag) {
	NLStructuralElement :: computeFirstPKStressVector(answer, gp, tStep);
      } else {
	FloatArray u, uN, uI;
	this->computeVectorOf({D_u, D_v, D_w}, VM_Total, tStep, u); // solution vector
	this->computeVectorOf({D_u, D_v, D_w}, VM_Total, tStep->givePreviousStep(), uN); // solution vector
	this->computeVectorOf({D_u, D_v, D_w}, VM_Incremental, tStep, uI); // solution vector
	double J0_J = computeFbarDeformationGradientVector(vF, gp, tStep, this, VM_Total);
	cs->giveFirstPKStresses(answer, gp, vF, tStep);
	answer.times(pow(J0_J,-2./3.));
      }
}


  

  
void
Structural2DElement :: giveMaterialOrientationAt(FloatArray &x, FloatArray &y, const FloatArray &lcoords)
{
    if ( this->elemLocalCS.isNotEmpty() ) { // User specified orientation
        x = {
            elemLocalCS.at(1, 1), elemLocalCS.at(2, 1)
        };
        y = {
            -x(1), x(0)
        };
    } else {
        FloatMatrix jac;
        this->giveInterpolation()->giveJacobianMatrixAt( jac, lcoords, * this->giveCellGeometryWrapper() );
        x.beColumnOf(jac, 1); // This is {dx/dxi, dy/dxi, dz/dxi}
        x.normalize();
        y = {
            -x(1), x(0)
        };
    }
}


double
Structural2DElement :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}




// Edge support

void
Structural2DElement :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */
    IntArray eNodes;
    static_cast< FEInterpolation2d * >( this->giveInterpolation() )->computeLocalEdgeMapping(eNodes,  iEdge);

    answer.resize(eNodes.giveSize() * 2);
    for ( int i = 1; i <= eNodes.giveSize(); i++ ) {
        answer.at(i * 2 - 1) = eNodes.at(i) * 2 - 1;
        answer.at(i * 2) = eNodes.at(i) * 2;
    }
}


double
Structural2DElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    /* Returns the line element ds associated with the given gp on the specific edge.
     * Note: The name is misleading since there is no volume to speak of in this case.
     * The returned value is used for integration of a line integral (external forces).
     */
    double detJ = static_cast< FEInterpolation2d * >( this->giveInterpolation() )->
                  edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );
    return detJ * gp->giveWeight();
}


int
Structural2DElement :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    FloatArray normal(2);

    static_cast< FEInterpolation2d * >( this->giveInterpolation() )->
    edgeEvalNormal( normal, iEdge, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );

    answer.resize(2, 2);
    answer.zero();
    answer.at(1, 1) = normal.at(2);
    answer.at(1, 2) = normal.at(1);
    answer.at(2, 1) = -normal.at(1);
    answer.at(2, 2) = normal.at(2);

    return 1;
}



// Edge support
void
Structural2DElement :: computeEdgeNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    /* Returns the [ 3 x (nno*3) ] shape function matrix {N} of the receiver, 
     * evaluated at the given gp.
     * {u} = {N}*{a} gives the displacements at the integration point.
     */ 
          
    // Evaluate the shape functions at the position of the gp. 
    FloatArray N;
    static_cast< FEInterpolation2d* > ( this->giveInterpolation() )->
        edgeEvalN( N, iSurf, sgp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );  
    answer.beNMatrixOf(N, 2);
}



// support for pressure follower load interface
void
Structural2DElement ::  surfaceEvaldNdxi(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{

  FloatMatrix dN;
  FEInterpolation2d *interp2d = static_cast<FEInterpolation2d*> (this->giveInterpolation());
  interp2d->edgeEvaldNdxi(dN, iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));
  answer.resize(2, 6);
  answer.zero();
  answer.at(1, 1) = dN.at(1,1);
  answer.at(1, 3) = dN.at(2,1);
  answer.at(1, 5) = dN.at(3,1);
  answer.at(2, 2) = dN.at(1,1);
  answer.at(2, 4) = dN.at(2,1);
  answer.at(2, 6) = dN.at(3,1);
  

}


void
Structural2DElement ::  surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdksi, int iSurf, GaussPoint *gp, TimeStep *tStep)
{

  answer.resize(2);
  IntArray bNodes;
  FloatArray gcoords, lcoords, vU, x;
  FloatMatrix dNdxi;

  lcoords = gp->giveNaturalCoordinates();
  this->giveBoundaryEdgeNodes(bNodes, iSurf);

  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);
  double nNodes = bNodes.giveSize();
  x.resize(2*nNodes);
  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    // compute actual node positions for Total Lagrangean formulation
    this->computeBoundaryVectorOf(bNodes, {D_u, D_v}, VM_Total, tStep, vU); // solution vector    
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i-1) = node->giveCoordinate(1) + vU.at( (i-1) * 2 + 1);
      x.at(2*i) = node->giveCoordinate(2) + vU.at( (i-1) * 2 + 2);
    }
  } else {
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i-1) = node->giveCoordinate(1);
      x.at(2*i) = node->giveCoordinate(2);
    }
    
  }
  
  FloatArray dx;
  FloatMatrix e3;
  e3 = {{0,-1},{1,0}};
  dx.beProductOf(dNdxi,x);
  answer.beProductOf(e3,dx);

}


void
Structural2DElement ::  surfaceEvalNormalAt(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdksi, int iSurf, GaussPoint *gp, TimeStep *tStep)
{

  answer.resize(2);
  IntArray bNodes;
  FloatArray gcoords, lcoords, vU, x;
  FloatMatrix dNdxi;

  lcoords = gp->giveNaturalCoordinates();
  this->giveBoundaryEdgeNodes(bNodes, iSurf);

  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);
  double nNodes = bNodes.giveSize();
  x.resize(2*nNodes);
  // compute actual node positions for Total Lagrangean formulation
  this->computeBoundaryVectorOf(bNodes, {D_u, D_v}, VM_Total, tStep, vU); // solution vector    
  for(int i = 1; i <= nNodes; i++) {
    Node *node = this->giveNode(bNodes.at(i));
    x.at(2*i-1) = node->giveCoordinate(1);
    x.at(2*i) = node->giveCoordinate(2);
  }
    
  FloatArray dx;
  FloatMatrix e3;
  e3 = {{0,-1},{1,0}};
  dx.beProductOf(dNdxi,x);
  answer.beProductOf(e3,dx);

}



void
Structural2DElement ::  surfaceEvalDeformedNormalAt_fromU(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdksi, int iSurf, GaussPoint *gp, TimeStep *tStep, const FloatArray &vU)
{

  answer.resize(2);
  IntArray bNodes;
  FloatArray gcoords, lcoords, x;
  FloatMatrix dNdxi;

  lcoords = gp->giveNaturalCoordinates();
  this->giveBoundaryEdgeNodes(bNodes, iSurf);

  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);
  double nNodes = bNodes.giveSize();
  x.resize(2*nNodes);
  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    // compute actual node positions for Total Lagrangean formulation
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i-1) = node->giveCoordinate(1) + vU.at( (i-1) * 2 + 1);
      x.at(2*i) = node->giveCoordinate(2) + vU.at( (i-1) * 2 + 2);
    }
  } else {
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i-1) = node->giveCoordinate(1);
      x.at(2*i) = node->giveCoordinate(2);
    }
    
  }
  
  FloatArray dx;
  FloatMatrix e3;
  e3 = {{0,-1},{1,0}};
  dx.beProductOf(dNdxi,x);
  answer.beProductOf(e3,dx);

}



void
Structural2DElement ::  surfaceEvalNumericalStiffMatrixAt(FloatMatrix &answer, FloatMatrix &dNdx,FloatArray &dxdeta, FloatArray &dxdxi, int iSurf, GaussPoint *gp, TimeStep *tStep)
{

  double r;
  IntArray bNodes;
  FloatArray lcoords, vU, dx, x;  
  FloatMatrix dNdxi, N, e3;

  e3 = {{0,-1},{1,0}};
  
  lcoords = gp->giveNaturalCoordinates();
  this->giveBoundaryEdgeNodes (bNodes, iSurf);
  double nNodes = bNodes.giveSize();
  x.resize(2*nNodes);
  this->surfaceEvalNmatrixAt(N, iSurf, gp);
  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);
  this->computeBoundaryVectorOf(bNodes, {D_u, D_v}, VM_Total, tStep, vU); // solution vector

  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    for(int i = 1; i <= 3; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i-1) = node->giveCoordinate(1) + vU.at( (i-1) * 2 + 1);
      x.at(2*i) = node->giveCoordinate(2) + vU.at( (i-1) * 2 + 2);
    }
  } else {
    for(int i = 1; i <= 2; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i -1) = node->giveCoordinate(1);
      x.at(2*i) = node->giveCoordinate(2);
      // Evaluate radius after deformation
      r += x.at(2*i -1) * N.at(1, 2 * i - 1);
    }
  }

  
  FloatArray n;
  dx.beProductOf(dNdxi,x);
  n.beProductOf(e3, dx);

  answer.resize(6,6);
  double pert = 1.e-6;
  FloatArray v, xp(x), vUp(vU), dn;
  FloatMatrix Ki,  Kn(2,6), Ktest(6,6);

  this->surfaceEvalDeformedNormalAt(n, dxdeta, dxdxi, iSurf, gp, tStep);

  
  for (int i = 1; i<=6; i++) {
    
    xp.at(i) += pert;      
    dx.beProductOf(dNdxi,xp);
    v.beProductOf(e3, dx);
    v.subtract(n);
    v.times(1/pert);
    Ki.beTProductOf(N, v);
    //
    vUp.at(i) += pert;
    this->surfaceEvalDeformedNormalAt_fromU(dn, dxdeta, dxdxi, iSurf, gp, tStep, vUp);
    for(int j = 1; j <= 2; j++) {
      Kn.at(j,i) = (dn.at(j) - n.at(j)) / pert;
    }
      
    for (int k = 1; k <= 6; k++) {
      answer.at(k,i) = Ki.at(k,1);
    }
    xp = x;
    vUp = vU;
    }
  dNdx = Kn;
  Ktest.beTProductOf(N, Kn);



  

  /*
  for (int i = 1; i<=4; i++) {
    dx.beProductOf(dNdxi,x);
    v.beProductOf(e3, dx);
    double rp = r+ pert;
    v.times(rp);
    v.subtract(n);
    v.times(1/pert);
    Ki.beTProductOf(N, v);
    for (int k = 1; k <= 4; k++) {
      answer.at(k,i) += Ki.at(k,1);
    }
    xp = x;
  }
  */
  
}




// Plane stress

PlaneStressElement :: PlaneStressElement(int n, Domain *aDomain) :
    Structural2DElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{}


void
PlaneStressElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, int lowerIndx, int upperIndx)
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx;
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper(tStep) );

    answer.resize(3, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = dNdx.at(i, 1);
        answer.at(2, i * 2 - 0) = dNdx.at(i, 2);

        answer.at(3, 2 * i - 1) = dNdx.at(i, 2);
        answer.at(3, 2 * i - 0) = dNdx.at(i, 1);
    }
}


void
PlaneStressElement :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, double alpha)
{
    // Returns the [ 4 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    /// @todo not checked if correct

    FloatMatrix dNdx;
    this->giveInterpolation()->evaldNdx( dNdx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper(tStep, alpha) );

    answer.resize(4, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 2 * i - 1) = dNdx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dNdx.at(i, 2);     // dv/dy -2
        answer.at(3, 2 * i - 1) = dNdx.at(i, 2);     // du/dy -6
        answer.at(4, 2 * i - 0) = dNdx.at(i, 1);     // dv/dx -9
    }
}


void
PlaneStressElement :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    if ( this->matRotation ) {
        ///@todo This won't work properly with "useUpdatedGpRecord" (!)
        FloatArray x, y;
        FloatArray rotStrain, s;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        // Transform to material c.s.
        rotStrain = {
            e(0) * x(0) * x(0) + e(2) * x(0) * x(1) + e(1) * x(1) * x(1),
            e(0) * y(0) * y(0) + e(2) * y(0) * y(1) + e(1) * y(1) * y(1),
            2 * e(0) * x(0) * y(0) + 2 * e(1) * x(1) * y(1) + e(2) * ( x(1) * y(0) + x(0) * y(1) )
        };

        this->giveStructuralCrossSection()->giveRealStress_PlaneStress(s, gp, rotStrain, tStep);

        answer = {
            s(0) * x(0) * x(0) + 2 * s(2) * x(0) * y(0) + s(1) * y(0) * y(0),
            s(0) * x(1) * x(1) + 2 * s(2) * x(1) * y(1) + s(1) * y(1) * y(1),
            s(1) * y(0) * y(1) + s(0) * x(0) * x(1) + s(2) * ( x(1) * y(0) + x(0) * y(1) )
        };
    } else {
        this->giveStructuralCrossSection()->giveRealStress_PlaneStress(answer, gp, e, tStep);
    }
}


void
PlaneStressElement :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_PlaneStress(answer, rMode, gp, tStep);
    if ( this->matRotation ) {
        FloatArray x, y;
        FloatMatrix Q;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );

        Q = {
            { x(0) * x(0), x(1) * x(1), x(0) * x(1) },
            { y(0) * y(0), y(1) * y(1), y(0) * y(1) },
            { 2 * x(0) * y(0), 2 * x(1) * y(1), x(1) * y(0) + x(0) * y(1) }
        };
        answer.rotatedWith(Q, 't');
    }
}





// Plane strain


PlaneStrainElement :: PlaneStrainElement(int n, Domain *aDomain) :
    Structural2DElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{}


void
PlaneStrainElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, int lowerIndx, int upperIndx)
// Returns the [ 4 x (nno*2) ] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx;
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper(tStep) );


    answer.resize(4, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = dNdx.at(i, 1);
        answer.at(2, i * 2 - 0) = dNdx.at(i, 2);

        answer.at(4, 2 * i - 1) = dNdx.at(i, 2);
        answer.at(4, 2 * i - 0) = dNdx.at(i, 1);
    }
}


void
PlaneStrainElement :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, double alpha)
{
    // Returns the [ 5 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    /// @todo not checked if correct

    FloatMatrix dNdx;
    this->giveInterpolation()->evaldNdx( dNdx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper(tStep, alpha) );

    answer.resize(5, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 2 * i - 1) = dNdx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dNdx.at(i, 2);     // dv/dy -2
        answer.at(4, 2 * i - 1) = dNdx.at(i, 2);     // du/dy -6
        answer.at(5, 2 * i - 0) = dNdx.at(i, 1);     // dv/dx -9
    }
}


void
PlaneStrainElement :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    if ( this->matRotation ) {
        ///@todo This won't work properly with "useUpdatedGpRecord" (!)
        FloatArray x, y;
        FloatArray rotStrain, s;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        // Transform to material c.s.
        rotStrain = {
            e(0) * x(0) * x(0) + e(3) * x(0) * x(1) + e(1) * x(1) * x(1),
            e(0) * y(0) * y(0) + e(3) * y(0) * y(1) + e(1) * y(1) * y(1),
            e(2),
            2 * e(0) * x(0) * y(0) + 2 * e(1) * x(1) * y(1) + e(3) * ( x(1) * y(0) + x(0) * y(1) )
        };
        this->giveStructuralCrossSection()->giveRealStress_PlaneStrain(s, gp, rotStrain, tStep);
        answer = {
            s(0) * x(0) * x(0) + 2 * s(3) * x(0) * y(0) + s(1) * y(0) * y(0),
            s(0) * x(1) * x(1) + 2 * s(3) * x(1) * y(1) + s(1) * y(1) * y(1),
            s(2),
            y(1) * ( s(3) * x(0) + s(1) * y(0) ) + x(1) * ( s(0) * x(0) + s(3) * y(0) )
        };
    } else {
        this->giveStructuralCrossSection()->giveRealStress_PlaneStrain(answer, gp, e, tStep);
    }
}


void
PlaneStrainElement :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_PlaneStrain(answer, rMode, gp, tStep);
    if ( this->matRotation ) {
        FloatArray x, y;
        FloatMatrix Q;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        Q = {
            { x(0) * x(0), x(1) * x(1), 0, x(0) * x(1) },
            { y(0) * y(0), y(1) * y(1), 0, y(0) * y(1) },
            { 0, 0, 1, 0 },
            { 2 * x(0) * y(0), 2 * x(1) * y(1), 0, x(1) * y(0) + x(0) * y(1) }
        };

        answer.rotatedWith(Q, 't');
    }
}



// Axisymmetry

AxisymElement :: AxisymElement(int n, Domain *aDomain) :
    Structural2DElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    //nlGeometry = 0; // Geometrical nonlinearities disabled as default
}


double
AxisymElement :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForAxisymmElements(normalToCrackPlane);
}


double
AxisymElement :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
  // note: radius is accounted by interpolation (of Fei2d*Axi type)
  double determinant = fabs( static_cast< FEInterpolation2d * >( this->giveInterpolation() )->
                             giveTransformationJacobian( gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() ) );

  double weight = gp->giveWeight();
  return determinant * weight;
  
}


void
AxisymElement :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, int li, int ui)
//
// Returns the [ 6 x (nno*2) ] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
    FEInterpolation *interp = this->giveInterpolation();

    FloatArray N;
    interp->evalN( N, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper(tStep) );
    double r = 0.0;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        double x = this->giveNode(i)->giveCoordinate(1);
        r += x * N.at(i);
    }

    FloatMatrix dNdx;
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );
    answer.resize(6, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = dNdx.at(i, 1);
        answer.at(2, i * 2 - 0) = dNdx.at(i, 2);
        answer.at(3, i * 2 - 1) = N.at(i) / r;
        answer.at(6, 2 * i - 1) = dNdx.at(i, 2);
        answer.at(6, 2 * i - 0) = dNdx.at(i, 1);
    }
}


void
AxisymElement :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, double alpha)
// Returns the [ 9 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz = u/r, 0, 0, du/dy,  0, 0, dv/dx
///@todo not checked if correct, is dw/dz = u/r for large deformations? /JB
{
    FloatArray n;
    FloatMatrix dnx;
    FEInterpolation2d *interp = static_cast< FEInterpolation2d * >( this->giveInterpolation() );

    interp->evalN( n, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper(tStep, alpha) );
    interp->evaldNdx( dnx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper(tStep, alpha) );


    int nRows = dnx.giveNumberOfRows();
    answer.resize(9, nRows * 2);
    answer.zero();

    double r = 0., x;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        x =  this->giveCellGeometryWrapper(tStep, alpha)->giveVertexCoordinates(i)->at(1);
        r += x * n.at(i);
    }


    // mode is _3dMat !!!!!! answer.at(4,*), answer.at(5,*), answer.at(7,*), and answer.at(8,*) is zero
    for ( int i = 1; i <= nRows; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);     // du/dx
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);     // dv/dy
        answer.at(6, 2 * i - 1) = dnx.at(i, 2);     // du/dy
        answer.at(9, 2 * i - 0) = dnx.at(i, 1);     // dv/dx
    }

    for ( int i = 0; i < this->giveNumberOfDofManagers(); i++ ) {
        answer.at(3, 2 * i + 1) = n.at(i + 1) / r;
    }
}


double
AxisymElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    FloatArray c(2);
    double result = static_cast< FEInterpolation2d * >( this->giveInterpolation() )->
                    edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );


    return result * gp->giveWeight();
    // note: radius taken into account by Fei*Axi interpolation
}


void
AxisymElement :: computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 6) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}

void
AxisymElement :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    if ( this->matRotation ) {
        ///@todo This won't work properly with "useUpdatedGpRecord" (!)
        FloatArray x, y;
        FloatArray rotStrain, s;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        // Transform to material c.s.
        rotStrain = {
            e(0) * x(0) * x(0) + e(5) * x(0) * x(1) + e(1) * x(1) * x(1),
            e(0) * y(0) * y(0) + e(5) * y(0) * y(1) + e(1) * y(1) * y(1),
            e(2),
            e(4) * y(0) + e(3) * y(1),
            e(4) * x(0) + e(3) * x(1),
            2 * e(0) * x(0) * y(0) + 2 * e(1) * x(1) * y(1) + e(5) * ( x(1) * y(0) + x(0) * y(1) )
        };
        this->giveStructuralCrossSection()->giveRealStress_3d(s, gp, rotStrain, tStep);
        answer = {
            s(0) * x(0) * x(0) + 2 * s(5) * x(0) * y(0) + s(1) * y(0) * y(0),
            s(0) * x(1) * x(1) + 2 * s(5) * x(1) * y(1) + s(1) * y(1) * y(1),
            s(2),
            s(4) * x(1) + s(3) * y(1),
            s(4) * x(0) + s(3) * y(0),
            y(1) * ( s(5) * x(0) + s(1) * y(0) ) + x(1) * ( s(0) * x(0) + s(5) * y(0) )
        };
    } else {
        this->giveStructuralCrossSection()->giveRealStress_3d(answer, gp, e, tStep);
    }
}

void
AxisymElement :: computeConstitutiveMatrixAt(FloatMatrix &answer,
                                             MatResponseMode rMode, GaussPoint *gp,
                                             TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_3d(answer, rMode, gp, tStep);
    if ( this->matRotation ) {
        FloatArray x, y;
        FloatMatrix Q;

        this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        Q = {
            { x(0) * x(0), x(1) * x(1), 0, 0, 0, x(0) * x(1) },
            { y(0) * y(0), y(1) * y(1), 0, 0, 0, y(0) * y(1) },
            { 0, 0, 1, 0, 0, 0 },
            { 0, 0, 0, y(1), y(0), 0 },
            { 0, 0, 0, x(1), x(0), 0 },
            { 2 * x(0) * y(0), 2 * x(1) * y(1), 0, 0, 0, x(1) * y(0) + x(0) * y(1) }
        };

        answer.rotatedWith(Q, 't');
    }
}
} // end namespace oofem
