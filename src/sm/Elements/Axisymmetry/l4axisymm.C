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

#include "../sm/Elements/Axisymmetry/l4axisymm.h"
#include "fei2dquadlin.h"
#include "node.h"
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

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(L4Axisymm);

FEI2dQuadLinAxi L4Axisymm :: interpolation(1, 2);

L4Axisymm :: L4Axisymm(int n, Domain *aDomain) :
    AxisymElement(n, aDomain), ZZNodalRecoveryModelInterface(this), SpatialLocalizerInterface(this)
{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 4;
    numberOfFiAndShGaussPoints = 1;
}


L4Axisymm :: ~L4Axisymm()
{ }


FEInterpolation *
L4Axisymm :: giveInterpolation() const { return & interpolation; }


Interface *
L4Axisymm :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    }
    
    return NULL;
}


IRResultType
L4Axisymm :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 4;
    IRResultType result = Structural2DElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }


    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 9 ) ||
           ( numberOfGaussPoints == 16 ) ) ) {
        numberOfGaussPoints = 4;
    }

    numberOfFiAndShGaussPoints = 4;

    return IRRT_OK;
}



void
L4Axisymm :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep, int li, int ui)
{
    // Returns the [ 6 x (nno*2) ] strain-displacement matrix {B} of the receiver,
    // evaluated at gp. Uses reduced integration.
    // (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
    // r = ( u1,v1,u2,v2,u3,v3,u4,v4)

    FloatArray N, NRed, redCoord;
    if ( numberOfFiAndShGaussPoints == 1 ) { // Reduced integration
        redCoord  = {0.0, 0.0}; // eval in centroid
    } else {
        redCoord = gp->giveNaturalCoordinates();
    }


    FEInterpolation *interp = this->giveInterpolation();
        

    interp->evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    interp->evalN( NRed, redCoord, FEIElementGeometryWrapper(this) );
    
    // Evaluate radius at center
    double r = 0.0;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        double x =  this->giveCellGeometryWrapper(tStep)->giveVertexCoordinates(i)->at(1);
        r += x * NRed.at(i);
    } 
    
    FloatMatrix dNdx, dNdxRed;
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper(tStep));
    interp->evaldNdx( dNdxRed, redCoord, *this->giveCellGeometryWrapper(tStep));
    answer.resize(6, dNdx.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = dNdx.at(i, 1);
        answer.at(2, i * 2 - 0) = dNdx.at(i, 2);
        answer.at(3, i * 2 - 1) = NRed.at(i) / r;
        answer.at(6, 2 * i - 1) = dNdxRed.at(i, 2);
        answer.at(6, 2 * i - 0) = dNdxRed.at(i, 1);
    }
}




  
  
  
void
L4Axisymm :: surfaceEvalNmatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{
  double ksi, n1, n2;
  
  ksi = gp->giveNaturalCoordinates().at(1);
  n1  = ( 1. - ksi ) * 0.5;
  n2  = ( 1. + ksi ) * 0.5;
  
  answer.resize(2, 4);
  answer.zero();
  answer.at(1, 1) = n1;
  answer.at(1, 3) = n2;
  answer.at(2, 2) = n1;
  answer.at(2, 4) = n2;
 

}

void
L4Axisymm :: surfaceEvaldNdxi(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{

  double dn1, dn2;
  
  dn1  =  - 0.5;
  dn2  =    0.5;
  
  answer.resize(2, 4);
  answer.zero();
  answer.at(1, 1) = dn1;
  answer.at(1, 3) = dn2;
  answer.at(2, 2) = dn1;
  answer.at(2, 4) = dn2;
  

}
  

IntegrationRule*
L4Axisymm :: surfaceGiveIntegrationRule(int order, int iSurf)
{
  return this->giveBoundaryEdgeIntegrationRule(order, iSurf);
}

void
L4Axisymm :: surfaceEvalNormalDerivative(FloatMatrix &answer, int iSurf, GaussPoint *gp, TimeStep *tStep)
{
  double r = 0.0;
  IntArray bNodes;
  FloatArray vU, x, dx;
  FloatMatrix e3, N, dNdxi;
 
  e3 = {{0,-1},{1,0}};
  this->giveBoundaryEdgeNodes (bNodes, iSurf);
  double nNodes = bNodes.giveSize();
  x.resize(2*nNodes);
  this->surfaceEvalNmatrixAt(N, iSurf, gp);
  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);


  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    // compute actual node positions for Total Lagrangean formulation
    this->computeBoundaryVectorOf(bNodes, {D_u, D_v}, VM_Total, tStep, vU); // solution vector    
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i-1) = node->giveCoordinate(1) + vU.at( (i-1) * 2 + 1);
      x.at(2*i) = node->giveCoordinate(2) + vU.at( (i-1) * 2 + 2);
      // Evaluate radius after deformation
      r += x.at(2*i-1) * N.at(1, 2 * i - 1);
    }
  } else {
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i -1) = node->giveCoordinate(1);
      x.at(2*i) = node->giveCoordinate(2);
      // Evaluate radius after deformation
      r += x.at(2*i-1) * N.at(1, 2 * i - 1);
    }
  }

  
  dx.beProductOf(dNdxi,x);
  FloatArray I(2);
  I.at(1) = 1.;

  FloatArray NI;
  FloatMatrix a1,a2, dxNI;

  
  
  a1.beProductOf(e3, dNdxi);
  a1.times(r);

  NI.beTProductOf(N,I);
  dxNI.beDyadicProductOf(dx,NI);
  a2.beProductOf(e3,dxNI);

  answer = a1;
  //answer.subtract(a2);
  answer.add(a2);
  /*answer = a2;
  answer.times(-1);
  */
}
  

  
void
L4Axisymm ::  surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdxi, int iSurf, GaussPoint *gp, TimeStep *tStep)
{

  IntArray bNodes;
  FloatArray vU, vdx, x, dx;
  FloatMatrix e3, N, dNdxi;
 
  e3 = {{0,-1},{1,0}};
  this->giveBoundaryEdgeNodes (bNodes, iSurf);
  double nNodes = bNodes.giveSize();
  x.resize(2*nNodes);
  this->surfaceEvalNmatrixAt(N, iSurf, gp);
  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);

  double r = 0.0;

  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    // compute actual node positions for Total Lagrangean formulation
    this->computeBoundaryVectorOf(bNodes, {D_u, D_v}, VM_Total, tStep, vU); // solution vector    
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i-1) = node->giveCoordinate(1) + vU.at( (i-1) * 2 + 1);
      x.at(2*i) = node->giveCoordinate(2) + vU.at( (i-1) * 2 + 2);
      // Evaluate radius after deformation
      r += x.at(2*i -1) * N.at(1, 2 * i - 1);
    }
  } else {
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i -1) = node->giveCoordinate(1);
      x.at(2*i) = node->giveCoordinate(2);
      // Evaluate radius after deformation
      r += x.at(2*i -1) * N.at(1, 2 * i - 1);
    }
  }

  
  dx.beProductOf(dNdxi,x);
  answer.beProductOf(e3, dx);
  answer.times(r);
  
}
  

void
L4Axisymm ::  surfaceEvalNormalAt(FloatArray &answer, FloatArray &dxdeta, FloatArray &dxdxi, int iSurf, GaussPoint *gp, TimeStep *tStep)
{

  IntArray bNodes;
  FloatArray vU, vdx, x, dx;
  FloatMatrix e3, N, dNdxi;
 
  e3 = {{0,-1},{1,0}};
  this->giveBoundaryEdgeNodes (bNodes, iSurf);
  double nNodes = bNodes.giveSize();
  x.resize(2*nNodes);
  this->surfaceEvalNmatrixAt(N, iSurf, gp);
  this->surfaceEvaldNdxi(dNdxi, iSurf, gp);

  double r = 0.0;

  if(this->domain->giveEngngModel()->giveFormulation() != AL) {
    // compute actual node positions for Total Lagrangean formulation
    this->computeBoundaryVectorOf(bNodes, {D_u, D_v}, VM_Total, tStep, vU); // solution vector    
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i-1) = node->giveCoordinate(1);
      x.at(2*i) = node->giveCoordinate(2);
      // Evaluate radius deformation
      r += x.at(2*i -1) * N.at(1, 2 * i - 1);
    }
  } else {
    for(int i = 1; i <= nNodes; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i -1) = node->giveCoordinate(1);
      x.at(2*i) = node->giveCoordinate(2);
      // Evaluate radius after deformation
      r += x.at(2*i -1) * N.at(1, 2 * i - 1);
    }
  }

  
  dx.beProductOf(dNdxi,x);
  answer.beProductOf(e3, dx);
  answer.times(r);
  
}

  

void
L4Axisymm ::  surfaceEvalNumericalStiffMatrixAt(FloatMatrix &answer, FloatArray &dxdeta, FloatArray &dxdxi, int iSurf, GaussPoint *gp, TimeStep *tStep)
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
    for(int i = 1; i <= 2; i++) {
      Node *node = this->giveNode(bNodes.at(i));
      x.at(2*i-1) = node->giveCoordinate(1) + vU.at( (i-1) * 2 + 1);
      x.at(2*i) = node->giveCoordinate(2) + vU.at( (i-1) * 2 + 2);
      // Evaluate radius after deformation
      r += x.at(2*i -1) * N.at(1, 2 * i - 1);
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
  n.times(r);

  answer.resize(4,4);
  double pert = 1.e-9;
  FloatArray v, xp(x);
  FloatMatrix Ki;


  /*for (int i = 1; i<=4; i++) {
    xp.at(i) += pert;      
    dx.beProductOf(dNdxi,xp);
    v.beProductOf(e3, dx);
    v.times(r);
    v.subtract(n);
    v.times(1/pert);
    Ki.beTProductOf(N, v);
    for (int k = 1; k <= 4; k++) {
      answer.at(k,i) = Ki.at(k,1);
    }
    xp = x;
    }*/


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



  
}
  










  

void
L4Axisymm :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    for ( int i = 1; i < 5; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
L4Axisymm :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i < 5; i++ ) {
        if ( pap == this->giveNode(i)->giveNumber() ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}

int
L4Axisymm :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
}


SPRPatchType
L4Axisymm :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}



#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void L4Axisymm :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 4 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void L4Axisymm :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 4 ];
    GraphicObj *go;
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}



void L4Axisymm :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 4 ];
    GraphicObj *tr;
    FloatArray v [ 4 ];
    double s [ 4 ], defScale;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 4; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], gc.giveIntVarType(), gc.giveIntVarMode(), i, tStep);
        }

        if ( result != 4 ) {
            return;
        }

        indx = gc.giveIntVarIndx();

        for ( i = 1; i <= 4; i++ ) {
            s [ i - 1 ] = v [ i - 1 ].at(indx);
        }

        if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
            for ( i = 0; i < 4; i++ ) {
                if ( gc.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    defScale = gc.getDefScale();
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                    p [ i ].z = 0.;
                } else {
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                    p [ i ].z = 0.;
                }
            }

            //EASValsSetColor(gc.getYieldPlotColor(ratio));
            gc.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);

            /*
             * } else if (gc.getScalarAlgo() == SA_ISO_LINE) {
             *
             * EASValsSetColor(context.getActiveCrackColor());
             * EASValsSetLineWidth(OOFEG_ISO_LINE_WIDTH);
             *
             * for (i=0; i< 4; i++) {
             * if (gc.getInternalVarsDefGeoFlag()) {
             * // use deformed geometry
             * defScale = gc.getDefScale();
             * p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,defScale);
             * p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,defScale);
             * p[i].z = 0.;
             *
             * } else {
             * p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
             * p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
             * p[i].z = 0.;
             * }
             * }
             *
             * // isoline implementation
             * oofeg_drawIsoLinesOnQuad (p, s);
             *
             */
        }
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        if ( numberOfGaussPoints != 4 ) {
            return;
        }

        IntArray ind(4);
        WCRec pp [ 9 ];

        for ( i = 0; i < 4; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                pp [ i ].z = 0.;
            } else {
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                pp [ i ].z = 0.;
            }
        }

        for ( i = 0; i < 3; i++ ) {
            pp [ i + 4 ].x = 0.5 * ( pp [ i ].x + pp [ i + 1 ].x );
            pp [ i + 4 ].y = 0.5 * ( pp [ i ].y + pp [ i + 1 ].y );
            pp [ i + 4 ].z = 0.5 * ( pp [ i ].z + pp [ i + 1 ].z );
        }

        pp [ 7 ].x = 0.5 * ( pp [ 3 ].x + pp [ 0 ].x );
        pp [ 7 ].y = 0.5 * ( pp [ 3 ].y + pp [ 0 ].y );
        pp [ 7 ].z = 0.5 * ( pp [ 3 ].z + pp [ 0 ].z );

        pp [ 8 ].x = 0.25 * ( pp [ 0 ].x + pp [ 1 ].x + pp [ 2 ].x + pp [ 3 ].x );
        pp [ 8 ].y = 0.25 * ( pp [ 0 ].y + pp [ 1 ].y + pp [ 2 ].y + pp [ 3 ].y );
        pp [ 8 ].z = 0.25 * ( pp [ 0 ].z + pp [ 1 ].z + pp [ 2 ].z + pp [ 3 ].z );

        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            const FloatArray &gpCoords = gp->giveNaturalCoordinates();
            if ( ( gpCoords.at(1) > 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 0;
                ind.at(2) = 4;
                ind.at(3) = 8;
                ind.at(4) = 7;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 4;
                ind.at(2) = 1;
                ind.at(3) = 5;
                ind.at(4) = 8;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) < 0. ) ) {
                ind.at(1) = 5;
                ind.at(2) = 2;
                ind.at(3) = 6;
                ind.at(4) = 8;
            } else {
                ind.at(1) = 6;
                ind.at(2) = 3;
                ind.at(3) = 7;
                ind.at(4) = 8;
            }

            if ( giveIPValue(v [ 0 ], gp, gc.giveIntVarType(), tStep) == 0 ) {
                return;
            }

            indx = gc.giveIntVarIndx();

            for ( i = 1; i <= 4; i++ ) {
                s [ i - 1 ] = v [ 0 ].at(indx);
            }

            for ( i = 0; i < 4; i++ ) {
                p [ i ].x = pp [ ind.at(i + 1) ].x;
                p [ i ].y = pp [ ind.at(i + 1) ].y;
                p [ i ].z = pp [ ind.at(i + 1) ].z;
            }

            gc.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}


#endif
} // end namespace oofem
