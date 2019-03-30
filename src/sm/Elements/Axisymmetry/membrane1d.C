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

#include "membrane1d.h"
#include "fei2dlinelin.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#include "mathfem.h"
#include "crosssection.h"
#include "classfactory.h"
#include "material.h"
#include "boundaryload.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/Materials/structuralms.h"

namespace oofem {
REGISTER_Element(Membrane1d);

FEI2dLineLin Membrane1d :: interp(1, 2);

Membrane1d :: Membrane1d(int n, Domain *aDomain) :
  NLStructuralElement(n, aDomain), PressureFollowerLoadElementInterface(this)
{
    numberOfDofMans = 2;
    nlGeometry = 1;
    pitch = 10;
}


Membrane1d :: ~Membrane1d()
{ }

void
Membrane1d :: postInitialize()
{
    // Element must be created before giveNumberOfNodes can be called
    StructuralElement :: postInitialize();
}


IRResultType
Membrane1d :: initializeFrom(InputRecord *ir)
{

    IRResultType result = NLStructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    initialStretch = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, initialStretch, _IFT_Membrane1d_InitialStretch);

    return IRRT_OK;
}




  

FEInterpolation *
Membrane1d :: giveInterpolation() const
{
    return & interp;
}

void 
Membrane1d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


double
Membrane1d :: computeLength()
// Returns the length of the receiver.
{

    double dx, dz;
    Node *nodeA, *nodeB;
    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dz      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        length  = sqrt(dx * dx + dz * dz);
    }

    return length;
}

void
Membrane1d :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{

    FloatArray N;
    interp.evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );


    answer.resize(2, 4);
    answer.zero();
    answer.at(1, 1) = N.at(1);
    answer.at(1, 3) = N.at(2);
    answer.at(2, 2) = N.at(1);
    answer.at(2, 4) = N.at(2);
}


void
Membrane1d :: computeNuMatrixAt(const FloatArray &iLocCoord, FloatArray &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{

    FloatArray N;
    interp.evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );

    answer.resize(4);
    answer.zero();
    answer.at(1) = N.at(1);
    answer.at(3) = N.at(2);
}


void
Membrane1d :: computeNwMatrixAt(const FloatArray &iLocCoord, FloatArray &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{

    FloatArray N;
    interp.evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );


    answer.resize(4);
    answer.zero();
    answer.at(2) = N.at(1);
    answer.at(4) = N.at(2);
}

void
Membrane1d :: computeBuMatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *tStep)
{

  FloatMatrix dNdx;
  interp.evaldNdxi( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
  double j  = interp.giveTransformationJacobian(gp->giveNaturalCoordinates(),FEIElementGeometryWrapper(this));
  dNdx.times(1./j);
  answer.resize(4);
  answer.zero();
  answer.at(1) = dNdx.at(1,1);
  answer.at(3) = dNdx.at(2,1);
}



void
Membrane1d :: computeBwMatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *tStep)
{

  FloatMatrix dNdx;
  interp.evaldNdxi( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
  double j  = interp.giveTransformationJacobian(gp->giveNaturalCoordinates(),FEIElementGeometryWrapper(this));
  dNdx.times(1./j);
  answer.resize(4);
  answer.zero();
  answer.at(2) = dNdx.at(1,1);
  answer.at(4) = dNdx.at(2,1);
}

void
Membrane1d :: computeBrMatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *tStep, bool linearized)
{
    // Returns the [ 6 x (nno*2) ] strain-displacement matrix {B} of the receiver,
    // evaluated at gp. Uses reduced integration
    // (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
    // r = ( u1,v1,u2,v2,u3,v3,u4,v4)



    FloatMatrix dNdx;
    interp.evaldNdxi( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    double j  = interp.giveTransformationJacobian(gp->giveNaturalCoordinates(),FEIElementGeometryWrapper(this));
    dNdx.times(1./j);
    FloatArray b1, b2;
    FloatMatrix A, A2;
    
    
    b1.resize(4);
    b1.zero();
    b1.at(1) = dNdx.at(1,1);
    b1.at(3) = dNdx.at(2,1);

    b2.resize(4);
    b2.zero();
    b2.at(2) = dNdx.at(1,1);
    b2.at(4) = dNdx.at(2,1);


    A2.beDyadicProductOf(b1,b1);
    A.beDyadicProductOf(b2,b2);
    A.add(A2);

    FloatArray d;
    this->computeVectorOf({D_u, D_v}, VM_Total, tStep, d); // solution vector

    
    if(d.giveSize()) {
      answer.beProductOf(A,d);
      if(!linearized) {
	answer.times(0.5);
      }
      answer.add(b1);
    } else {
      answer = b1;
    }

}

void
Membrane1d :: computeBtMatrixAt(GaussPoint *gp, FloatArray &answer, TimeStep *tStep, bool linearized)
{
    // Returns the [ 6 x (nno*2) ] strain-displacement matrix {B} of the receiver,
    // evaluated at gp. Uses reduced integration
    // (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
    // r = ( u1,v1,u2,v2,u3,v3,u4,v4)


    FloatArray N;
    interp.evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    FloatArray n1;
    this->computeNuMatrixAt(gp->giveNaturalCoordinates(), n1);
         
    FloatMatrix dNdx;
    interp.evaldNdxi( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    double j  = interp.giveTransformationJacobian(gp->giveNaturalCoordinates(),FEIElementGeometryWrapper(this));
    dNdx.times(1./j);
    // Evaluate radius
    double r = 0.0;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
      double x = this->giveNode(i)->giveCoordinate(1);
      r += x * N.at(i);
    } 

    FloatArray d;
    this->computeVectorOf({D_u, D_v}, VM_Total, tStep, d); // solution vector
    double dn1;

    if(d.giveSize()) {
      if(linearized) {
	dn1 = 1./r + n1.dotProduct(d)/r/r;
      } else {
	dn1 = 1./r + 0.5*n1.dotProduct(d)/r/r;
      }      
    } else {
      dn1 = 1./r;
    }

    answer = n1;
    answer.times(dn1);
}




void
Membrane1d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, bool linearized)
{
    // Returns the [ 6 x (nno*2) ] strain-displacement matrix {B} of the receiver,
    // evaluated at gp. Uses reduced integration
    // (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
    // r = ( u1,v1,u2,v2,u3,v3,u4,v4)
    answer.resize(2,4);
    FloatArray br,bt;
    computeBrMatrixAt(gp, br, tStep, linearized);
    computeBtMatrixAt(gp, bt, tStep, linearized);

    answer.addSubVectorRow(br, 1, 1);
    answer.addSubVectorRow(bt, 2, 1);


}



void
Membrane1d :: computeG1G2MatricesAt(FloatMatrix &G1, FloatMatrix &G2,GaussPoint *gp, TimeStep *tStep)
{
  
  FloatArray n1, b1, b2;
  FloatMatrix A1;

  this->computeNuMatrixAt(gp->giveNaturalCoordinates(), n1);
  this->computeBuMatrixAt(gp, b1, tStep);
  this->computeBwMatrixAt(gp, b2, tStep);

  A1.beDyadicProductOf(b1,b1);
  G1.beDyadicProductOf(b2,b2);
  G1.add(A1);

  // Evaluate radius
  FloatArray N;
  interp.evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
  double r = 0.0;
  for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
    double x = this->giveNode(i)->giveCoordinate(1);
    r += x * N.at(i);
  } 

  G2.beDyadicProductOf(n1,n1);
  G2.times(1./r/r);
}

void
Membrane1d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v
    };
}


double
Membrane1d :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
  FloatArray N;
  interp.evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
  double r = 0.0;
  for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
    double x  = this->giveNode(i)->giveCoordinate(1);
    r += x * N.at(i);
  }
  double weight  = gp->giveWeight();
  return  interp.giveTransformationJacobian(gp->giveNaturalCoordinates(),FEIElementGeometryWrapper(this)) * weight * r * this->giveStructuralCrossSection()->give(CS_Thickness, gp);
  // * this->giveCrossSection()->give(CS_Area, gp);
}

bool
Membrane1d :: giveRotationMatrix(FloatMatrix &answer)
{
  answer.clear();
  return false;  
}

void
Membrane1d :: computeStiffnessMatrix(FloatMatrix &answer,MatResponseMode rMode, TimeStep *tStep)
{

  answer.clear();


  double dV;
  FloatArray u, strain, S;
  FloatMatrix bE, bj, d, dbj;
  //  bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
  for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
    this->computeBmatrixAt(gp, bj, tStep, true);
    this->computeBmatrixAt(gp, bE, tStep, false);
    this->giveStructuralCrossSection()->giveStiffnessMatrix_AxisymMembrane1d(d, rMode, gp, tStep);
    this->computeVectorOf(VM_Total, tStep, u);
    strain.beProductOf(bE, u);
    // add influence of initial stress/stretch
    double l2 = initialStretch*initialStretch;
    strain.times(l2);
    FloatArray E0(2);
    E0.at(1) = (l2-1.)/2;
    E0.at(2) = (l2-1.)/2;
    strain.add(E0);
    /////////////////////////////////////////////////////////////////////////////////////////
    this->giveStructuralCrossSection()->giveRealStress_AxisymMembrane1d(S, gp, strain, tStep);
    dV = this->computeVolumeAround(gp);
    dbj.beProductOf(d, bj);

    FloatMatrix G1,G2;
    this->computeG1G2MatricesAt(G1,G2,gp, tStep);

    /*
    if ( matStiffSymmFlag ) {
      answer.plusProductSymmUpper(bj, dbj, dV);  

    } else {
    */
      answer.plusProductUnsym(bj, dbj, dV);    
      // }

    if(S.at(1) == 0 && S.at(2) == 0 && initialStretch == 1 ) {
      S.at(1) = 0.1;
      S.at(2) = 0.1;
    }
    answer.add(S.at(1)*dV,G1);
    answer.add(S.at(2)*dV,G2);
    

  }
  /*
  if ( matStiffSymmFlag ) {
    answer.symmetrized();
  }
  */

}

void
Membrane1d :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  FloatMatrix b, bE;
  FloatArray u, stress, strain;
  // This function can be quite costly to do inside the loops when one has many slave dofs.
  this->computeVectorOf(VM_Total, tStep, u);
  // subtract initial displacements, if defined
  if ( initialDisplacements ) {
    u.subtract(* initialDisplacements);
  }

  // zero answer will resize accordingly when adding first contribution
  answer.clear();

  for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
    this->computeBmatrixAt(gp, b, tStep, true);
    this->computeBmatrixAt(gp, bE, tStep,false);
    
    if ( !this->isActivated(tStep) ) {
      strain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
      strain.zero();
    }
    strain.beProductOf(bE, u);
    // add influence of initial stress/stretch
    double l2 = initialStretch*initialStretch;
    strain.times(l2);
    FloatArray E0(2);
    E0.at(1) = (l2-1.)/2.;
    E0.at(2) = (l2-1.)/2.;
    strain.add(E0);
    //    b.times(l2);
    /////////////////////////////////////////////////////////////////////////////////////////
    this->giveStructuralCrossSection()->giveRealStress_AxisymMembrane1d(stress, gp, strain, tStep);    

    // updates gp stress and strain record  acording to current
    // increment of displacement
    if ( stress.giveSize() == 0 ) {
      break;
    }
    // compute nodal representation of internal forces using f = B^T*Sigma dV
    double dV = this->computeVolumeAround(gp);
    answer.plusProduct(b, stress, dV);      
  }
  
  // if inactive update state, but no contribution to global system
  if ( !this->isActivated(tStep) ) {
    answer.zero();
    return;
  }
  
}




void
Membrane1d :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->computeGlobalCoordinates( answer, gp->giveNaturalCoordinates() );
}

void
Membrane1d :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    /*
     *
     * computes interpolation matrix for element edge.
     * we assemble locally this matrix for only nonzero
     * shape functions.
     * (for example only two nonzero shape functions for 2 dofs are
     * necessary for linear plane stress tringle edge).
     * These nonzero shape functions are then mapped to
     * global element functions.
     *
     * Using mapping technique will allow to assemble shape functions
     * without regarding particular side
     */

    this->computeNmatrixAt(gp->giveSubPatchCoordinates(), answer);
}


void
Membrane1d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge != 1 ) {
        OOFEM_ERROR("wrong edge number");
    }


    answer.resize(4);
    answer.at(1) = 1;
    answer.at(2) = 2;
    answer.at(3) = 3;
    answer.at(4) = 4;
}

double
Membrane1d ::   computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
      OOFEM_ERROR("wrong egde number");
    }
    FloatArray N;
    interp.evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    double weight  = gp->giveWeight();
    return  interp.giveTransformationJacobian(gp->giveNaturalCoordinates(),FEIElementGeometryWrapper(this)) * weight;
 
}




int
Membrane1d :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    //double dx,dy, length ;
    double sine, cosine;

    answer.resize(2, 2);
    answer.zero();

    sine           = sin( this->givePitch() );
    cosine         = cos(pitch);

    answer.at(1, 1) = cosine;
    answer.at(1, 2) = -sine;
    answer.at(2, 1) = sine;
    answer.at(2, 2) = cosine;

    return 1;
}

double Membrane1d :: givePitch()
// Returns the pitch of the receiver.
{
    double xA, xB, zA, zB;
    Node *nodeA, *nodeB;

    if ( pitch == 10. ) {             // 10. : dummy initialization value
        nodeA  = this->giveNode(1);
        nodeB  = this->giveNode(2);
        xA     = nodeA->giveCoordinate(1);
        xB     = nodeB->giveCoordinate(1);
        zA     = nodeA->giveCoordinate(2);
        zB     = nodeB->giveCoordinate(2);
        pitch  = atan2(zB - zA, xB - xA);
    }

    return pitch;
}




bool
Membrane1d :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver.
{
    double sine, cosine;

    answer.resize(4, 4);
    answer.zero();

    sine = sin( this->givePitch() );
    cosine  = cos(pitch);
    answer.at(1, 1) =  cosine;
    answer.at(1, 2) =  sine;
    answer.at(2, 1) = -sine;
    answer.at(2, 2) =  cosine;

    answer.at(3, 3) =  cosine;
    answer.at(3, 4) =  sine;
    answer.at(4, 3) = -sine;
    answer.at(4, 4) =  cosine;

    return true;
}



void
Membrane1d :: surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdksi, FloatArray &dxdeta, int iSurf, GaussPoint *gp, TimeStep *tStep)
 {

     answer.resize(2);
     const FloatArray &lcoords = gp->giveNaturalCoordinates();
     FEInterpolation *fei = this->giveInterpolation();

     FloatArray vN;
     fei->evalN( vN, lcoords, FEIElementGeometryWrapper(this) );
     // Evaluate radius
     double r = 0.0;
     for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
	  double x = this->giveNode(i)->giveCoordinate(1);
	  r += x * vN.at(i);
     } 

     
     FloatArray Bu,Bw,Nu;
     this->computeNuMatrixAt(lcoords, Nu);
     this->computeBuMatrixAt(gp, Bu, tStep);
     this->computeBwMatrixAt(gp, Bw, tStep);

     FloatArray d;
     this->computeVectorOf({D_u, D_v}, VM_Total, tStep, d); // solution vector
      
     double dw, u, du;       
     u = Nu.dotProduct(d);
     du = Bu.dotProduct(d);
     dw = Bw.dotProduct(d);
     
     answer.at(1) = -dw*(r+u);
     answer.at(2) = (r+u)*(1+du);

 }




void
Membrane1d :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep)
{
    answer.clear();

    if ( edge != 1 ) {
        OOFEM_ERROR("Membrane1d only has 1 edge (the midline) that supports loads. Attempted to apply load to edge %d", edge);
    }

    if ( type != ExternalForcesVector ) {
        return;
    }
    FloatArray coords, t;
    FloatMatrix N, T;
    
    FloatArray d;
    this->computeVectorOf({D_u, D_v}, VM_Total, tStep, d); // solution vector

    answer.clear();
    for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        this->computeNmatrixAt(lcoords, N);
	FloatArray vN;
	interp.evalN( vN, lcoords, FEIElementGeometryWrapper(this) );
	// Evaluate radius
	double r = 0.0;
	for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
	  double x = this->giveNode(i)->giveCoordinate(1);
	  r += x * vN.at(i);
	} 

        if ( load ) {
            this->computeGlobalCoordinates(coords, lcoords);
            load->computeValues(t, tStep, coords, { D_u, D_v}, mode);
        } else {
            load->computeValues(t, tStep, lcoords, { D_u, D_v}, mode);
        }

        if ( load->giveCoordSystMode() == Load :: CST_Global ) {
            if ( this->computeLoadGToLRotationMtrx(T) ) {
                t.rotatedWith(T, 'n');
            }
        }

	FloatArray Bu,Bw,Nu, Nw;
	this->computeNuMatrixAt(lcoords, Nu);
	this->computeNwMatrixAt(lcoords, Nw);
	this->computeBuMatrixAt(gp, Bu, tStep);
	this->computeBwMatrixAt(gp, Bw, tStep);
	
	double dw, u, du;       
	u = Nu.dotProduct(d);
	du = Bu.dotProduct(d);
	dw = Bw.dotProduct(d);
	
	FloatArray f(2);
	
	f.at(1) = -dw*(r+u)*t.at(2);
	f.at(2) = (r+u)*(1+du)*t.at(2);

	
        double dl = this->computeEdgeVolumeAround(gp,edge);

        answer.plusProduct(N, f, dl);
    }

    // Loads from sets expects global c.s.
    this->computeGtoLRotationMatrix(T);
    answer.rotatedWith(T, 't');
}



Interface *
Membrane1d :: giveInterface(InterfaceType interface)
{
    if ( interface == PressureFollowerLoadElementInterfaceType ) {
        return static_cast< PressureFollowerLoadElementInterface * >(this);
    }
    
    return NULL;
}


double
Membrane1d :: computeMidsurfaceVolume(TimeStep *tStep)
{

    FloatArray d;
    double V = 0;
    this->computeVectorOf({D_u, D_v}, VM_Total, tStep, d); // solution vector
    for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
      FloatMatrix N;
      const FloatArray &lcoords = gp->giveNaturalCoordinates();
      this->computeNmatrixAt(lcoords, N);
      FloatArray vN;
      this->giveInterpolation()->evalN( vN, lcoords, FEIElementGeometryWrapper(this) );
      // Evaluate radius
      double r = 0.0;
      for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
	double x = this->giveNode(i)->giveCoordinate(1);
	r += x * vN.at(i);
      } 

      //      double dV = this->computeVolumeAround(gp);

      FloatArray Bu,Nw,Nu;
      this->computeNuMatrixAt(lcoords, Nu);
      this->computeBuMatrixAt(gp, Bu, tStep);
      this->computeNwMatrixAt(lcoords, Nw);
      
      double w, u, du;       
      u = Nu.dotProduct(d);
      du = Bu.dotProduct(d);
      w = Nw.dotProduct(d);

      double dl = this->computeEdgeVolumeAround(gp,1);

      V += w*(r+u)*(1+du)*dl;
     

      

    }
    return V;

}






} // end namespace oofem
