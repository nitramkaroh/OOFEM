#include "../sm/Elements/Beams/nlbeam_sm.h"
#include "material.h"
#include "crosssection.h"
#include "node.h"

#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(NlBeam_SM);


NlBeam_SM :: NlBeam_SM(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    this->internalForces.resize(3);

}


void
NlBeam_SM :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_w, R_v
    };
}



double
NlBeam_SM :: computeLength()
// Returns the length of the receiver.
{
    double dx, dy;
    Node *nodeA, *nodeB;

    if ( beamLength == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        beamLength  = sqrt(dx * dx + dy * dy);
    }

    return beamLength;
}
  

 double
NlBeam_SM :: givePitch()
// Returns the pitch of the receiver.
{
    double xA, xB, yA, yB;
    Node *nodeA, *nodeB;

    if ( pitch == 10. ) {             // 10. : dummy initialization value
        nodeA  = this->giveNode(1);
        nodeB  = this->giveNode(2);
        xA     = nodeA->giveCoordinate(1);
        xB     = nodeB->giveCoordinate(1);
        yA     = nodeA->giveCoordinate(3);
        yB     = nodeB->giveCoordinate(3);
        pitch  = atan2(yB - yA, xB - xA);
    }

    return pitch;
}

 IRResultType
NlBeam_SM :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // first call parent
    StructuralElement :: initializeFrom(ir);

    // Numerical parameters
    // 1. number of segments for numerical integration along the beam, default value 100
    IR_GIVE_OPTIONAL_FIELD(ir, NIP, _IFT_NlBeam_SM_NIP);
    IR_GIVE_FIELD(ir, EA, _IFT_NlBeam_SM_EA);
    IR_GIVE_FIELD(ir, EI, _IFT_NlBeam_SM_EI);
    /*    // relative tolerance for iterations at the section level
    IR_GIVE_OPTIONAL_FIELD(ir, section_tol, _IFT_FbarElementExtensionInterface_fbarflag);
    // maximum number of iterations at the section level
    IR_GIVE_OPTIONAL_FIELD(ir, section_maxit, _IFT_FbarElementExtensionInterface_fbarflag);
    */
    
    // relative tolerance for iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_tol, _IFT_NlBeam_SM_Beam_Tolerance);
    // maximum number of iterations at the beam level
    IR_GIVE_OPTIONAL_FIELD(ir, beam_maxit, _IFT_NlBeam_SM_Beam_MaxIteration);
    IR_GIVE_OPTIONAL_FIELD(ir, nsubsteps_init, _IFT_NlBeam_SM_Beam_NumberMaxSubsteps);
    
    this->givePitch();
    this->computeLength();

    this->x.resize(NIP+1);
    this->u.resize(NIP+1);
    this->w.resize(NIP+1);
    this->phi.resize(NIP+1);

    this->vN.resize(NIP+1);
    this->vV.resize(NIP+1);
    this->vM.resize(NIP+1);

    
    return IRRT_OK;
}


 

 

double
NlBeam_SM :: computeMomentFromCurvature(double kappa)
{
   return EI*kappa;
}

double
NlBeam_SM :: computeDerMomentFromCurvature(double kappa)
{
  return EI;
}

double
NlBeam_SM :: computeCurvatureFromMoment(double M)
{
  return M/EI;
}


void
NlBeam_SM :: integrateAlongBeamAndGetJacobi(const FloatArray &fab, FloatArray &ub, FloatMatrix &jacobi)
{
  ub.resize(3);
  FloatArray du(3), dw(3), dphi(3), dM(3), dkappa(3), dphi_mid(3), dN_mid(3);
  double Xab = fab.at(1), Zab = fab.at(2), Mab = fab.at(3);
  this->vM.at(1) = -Mab + Xab*w.at(1) - Zab * u.at(1);
  this->vN.at(1) = -Xab * cos(phi.at(1)) + Zab * sin(phi.at(1));
  double dx = beamLength/NIP;
  for (int i=2; i <= NIP+1; i++) {
    this->x.at(i) = this->x.at(i-1) + dx;
    double M = -Mab + Xab*w.at(i-1) - Zab * (x.at(i-1) + u.at(i-1));
    dM = {w.at(i-1)+Xab*dw.at(1)-Zab*du.at(1), -x.at(i-1) - u.at(i-1)+Xab*dw.at(2)-Zab*du.at(2), -1.+Xab*dw.at(3) - Zab*du.at(3)};
    double kappa = computeCurvatureFromMoment(M);
    double dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa = dM;
    dkappa.times(1./dMdkappa);
    double phi_mid = phi.at(i-1) + kappa * dx/2.;
    dphi_mid = dphi;
    dphi_mid.add(dx/2 ,dkappa);
    double N_mid = -Xab * cos(phi_mid) + Zab * sin(phi_mid);
    vN.at(i) = N_mid;
    dN_mid = ( Xab * sin(phi_mid) + Zab * cos(phi_mid)) * dphi_mid;
    dN_mid.at(1) -= cos(phi_mid);
    dN_mid.at(2) += sin(phi_mid);
    u.at(i) = u.at(i-1) + dx * ( ( 1. + N_mid / EA ) * cos(phi_mid) - 1. );
    du.add((dx/EA) * cos(phi_mid),dN_mid);
    du.add(- dx * (1. + N_mid/EA) * sin(phi_mid), dphi_mid);
    w.at(i) = w.at(i-1) -  dx * ( 1. + N_mid / EA ) * sin(phi_mid);
    dw.add( -dx / EA * sin(phi_mid),dN_mid);
    dw.add( -dx * ( 1. + N_mid / EA) * cos(phi_mid), dphi_mid);
    M = -Mab + Xab * w.at(i) - Zab * ( x.at(i) + u.at(i) );
    vM.at(i) = M;
    dM = {w.at(i) + Xab*dw.at(1) - Zab * du.at(1), -x.at(i)-u.at(i) + Xab * dw.at(2) - Zab * du.at(2), -1. + Xab * dw.at(3) - Zab * du.at(3)};
    kappa = computeCurvatureFromMoment(M);
    dMdkappa = computeDerMomentFromCurvature(kappa);
    dkappa = dM;
    dkappa.times(1./dMdkappa);
    phi.at(i) = phi_mid + kappa * dx / 2.;
    dphi = dphi_mid;
    dphi.add(dx/2 ,dkappa);
  }
  ub.at(1) = u.at(NIP+1); ub.at(2) = w.at(NIP+1); ub.at(3) = phi.at(NIP+1);

  jacobi.copySubVectorRow(du, 1,1);
  jacobi.copySubVectorRow(dw, 2,1);
  jacobi.copySubVectorRow(dphi, 3,1);

}


/*
Find forces and moment at the left end that lead to given displacements and rotations
at the right end, assuming that the displacement and rotation at the left end are zero.
The solution is computed iteratively, using the Newton-Raphson method. 
Everything is done here in the local coordinate system aligned with the undeformed beam.
*/
bool
NlBeam_SM :: findLeftEndForcesLocal(FloatArray &ub_target, FloatArray &fab_loc)
{
  FloatArray res, dforces, ub_loc, f_init(fab_loc);
  FloatMatrix jacobi(3,3);
  int iter = 0;
   double tolerance = beam_tol;
  if(ub_target.computeNorm() > 1.e-10) {
    tolerance *= ub_target.computeNorm();
  }
  this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, jacobi);
  res = ub_target-ub_loc; 
  double error = res.computeNorm();
  
  while ((iter==0) || (error>tolerance && iter < beam_maxit) ){
    //|| (error != error)
    iter++;
    jacobi.solveForRhs(res, dforces);
    fab_loc.add(dforces);
    this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, jacobi);
    res = ub_target - ub_loc; 
    error = res.computeNorm();
  }

  if (iter >= beam_maxit || error != error) {
    fab_loc = f_init;
    return false;
  }
  return true;
}

/*
Auxiliary matrix for transformation of displacements (or in transposed form of forces).
It corresponds to rotation from global axes to the auxiliary system aligned with the deformed beam.
*/
void
NlBeam_SM :: construct_T(FloatMatrix &T, const double phia)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = cos(pitch-phia);
  T.at(1,2) = sin(pitch-phia);
  T.at(2,1) = - T.at(1,2);
  T.at(3,3) = 1.;
}




/*
Auxiliary matrix for transformation of stiffness.
It is obtained by differentiating T with respect to phia.
*/
void
  NlBeam_SM :: construct_Tprime(FloatMatrix &T, const double phia)
{
  T.resize(3,3);
  T.at(1,1) = T.at(2,2) = sin(pitch-phia);
  T.at(1,2) = -cos(pitch-phia);
  T.at(2,1) = -T.at(1,2);
}

void
NlBeam_SM :: construct_l(FloatArray &l, double phia)
{
  l.resize(3);
  l.at(1) = beamLength * (cos(phia)-1.);
  l.at(2) = beamLength * sin(phia);
}


void
NlBeam_SM :: construct_l(FloatArray &l, double phia, double L)
{
  l.resize(3);
  l.at(1) = L * (cos(phia)-1.);
  l.at(2) = L * sin(phia);
}


void
NlBeam_SM :: construct_lprime(FloatArray &l, const double phia)
{
  l.resize(3);
  l.at(1) = -beamLength * sin(phia);
  l.at(2) = beamLength * cos(phia);
}


/*
Find forces and moment at the left end that lead to given displacements and rotations
at both ends.
The input displacements and output forces are taken with respect to the global coordinate system.
Note that the transformation matrix T is affected by angle alpha that specifies the initial beam geometry. 
*/
void
NlBeam_SM :: findLeftEndForces(const FloatArray &u, FloatArray &fab)
{
  FloatArray ub_loc, l, fab_loc, fab_init(fab); 
  FloatMatrix T; 
  // compute displacements of the right end with respect to the auxiliary coordinate system
  construct_l(ub_loc, u.at(3));
  construct_T(T, u.at(3));

  FloatArray u_b, u_a, u_ba, temp;
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDifferenceOf(u_b, u_a);
  temp.beProductOf(T,u_ba);
  ub_loc.add(temp);
  // transform initial guess to local coordinates
  fab_loc.beProductOf(T, fab);
  // find end forces in the local coordinate syste
  bool success = findLeftEndForcesLocal(ub_loc, fab_loc);
  if (!success) { // multi-level substepping (note that fab_loc remained unchanged if the local iteration failed)
    FloatArray ub_loc_substep, ub_loc_intermediate(3);
    FloatMatrix jacobi(3,3);
    this->integrateAlongBeamAndGetJacobi(fab_loc, ub_loc_intermediate, jacobi);
    ub_loc_substep.beDifferenceOf(ub_loc, ub_loc_intermediate);
    int refinement_factor = 4;
    ub_loc_substep.times(1./refinement_factor);
    int isubstep = 0, maxnsubsteps = 10000;
    int nsubsteps = this->nsubsteps_init;
    while (isubstep < nsubsteps && nsubsteps <= maxnsubsteps){
      ub_loc_intermediate.add(ub_loc_substep);
      success = findLeftEndForcesLocal(ub_loc_intermediate, fab_loc);
      if (success){
	isubstep++;
      } else { // further refinement needed
	ub_loc_intermediate.subtract(ub_loc_substep);
	ub_loc_substep.times(1./refinement_factor);
	nsubsteps = isubstep + 4*(nsubsteps-isubstep);
      }
    }
    if (success == false){
      fab = fab_init;
      domain->giveEngngModel()->setAnalysisCrash(true);
      return;
    }
  }


    
  // transform local end forces to the global coordinate syste
  fab.beTProductOf(T, fab_loc);
}

/*
Find end forces and moments that lead to given displacements and rotations
at both ends.
The input displacements and output forces are taken with respect to the global coordinate system.
*/
 void
NlBeam_SM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  // solution vector
  answer.resize(6);
  FloatArray u, f_a(3);
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
  // only the first three entries of f are computed
  this->findLeftEndForces(u, this->internalForces);
  answer.at(1) = this->internalForces.at(1);
  answer.at(2) = this->internalForces.at(2);
  answer.at(3) = this->internalForces.at(3);  
  answer.at(4) = -answer.at(1);
  answer.at(5) = -answer.at(2);
  double c1 = beamLength*sin(pitch) + u.at(5) - u.at(2);
  double c2 = -beamLength*cos(pitch) - u.at(4) + u.at(1);
  answer.at(6) = c1 * answer.at(1) + c2 * answer.at(2) - answer.at(3);
}


/*
Evaluate the tangent stiffness matrix based on given displacements at both ends
and given end forces at the left end (must be provided, which means that typically findLeftEndForces is run first).
*/
void
NlBeam_SM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  FloatArray lprime, fab_loc, ub_loc, Tu;
  FloatMatrix T, Tprime, G, Ginv, TtGinv;
  answer.resize(6,6);
  // solution vector
  FloatArray u;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u); 
  // compute auxiliary matrices
  construct_T(T, u.at(3)); 
  construct_Tprime(Tprime, u.at(3)); 
  construct_lprime(lprime, u.at(3));
  // transform left-end forces to local coordinates
  FloatArray fab(3);
  fab = this->internalForces;
  //this->findLeftEndForces(u, fab);
  //fab_loc.beProductOf(T, this->internalForces);
  fab_loc.beProductOf(T, fab);

  // get Jacobi matrix in local coordinates (ub_loc is dummy, will not be used)
  integrateAlongBeamAndGetJacobi(fab_loc, ub_loc, G);

  Ginv.beInverseOf(G);
  // compute product Ttransposed*Ginverse
  TtGinv.beTProductOf(T,Ginv);
  // compute product Ttransposed*Ginverse*T and store it in the upper stiffness block 
  FloatMatrix aux;
  aux.beProductOf(TtGinv, T);
  answer.setSubMatrix(aux, 1,1);
  answer.times(-1);
  answer.setSubMatrix(aux, 1,4);
  
  // compute Tprime*(ub-ua)
  FloatArray u_b, u_a, u_ba;
  u_a.beSubArrayOf(u, {1,2,3});
  u_b.beSubArrayOf(u, {4,5,6});
  u_ba.beDifferenceOf(u_b, u_a);
  Tu.beProductOf(Tprime, u_ba);
  // compute additional stiffness terms in the third column (associated with phia)
  FloatArray col3_1,col3_2;
  Tu.add(lprime);
  col3_1.beProductOf(TtGinv,Tu);
  col3_2.beTProductOf(Tprime, fab_loc);
  answer.addSubVectorCol(col3_1, 1,3);
  answer.addSubVectorCol(col3_2, 1,3);
  // construct the sixth row (complete formula)
  double c1 = beamLength*sin(pitch) + u.at(5) - u.at(2);
  double c2 = -beamLength*cos(pitch) - u.at(4) + u.at(1);
  answer.at(6,1) =  fab.at(2);//this->internalForces.at(2);
  answer.at(6,2) = -fab.at(1);//-this->internalForces.at(1);
  answer.at(6,4) = -fab.at(2);//-this->internalForces.at(2);
  answer.at(6,5) =  fab.at(1);//this->internalForces.at(1);
 
  for(int j = 1; j <= 6; j++) {
    answer.at(4,j) = -answer.at(1,j);
    answer.at(5,j) = -answer.at(2,j);
    answer.at(6,j) += c1 * answer.at(1,j) + c2 * answer.at(2,j) - answer.at(3,j);
  }
  FloatMatrix test;
  //this->computeStiffnessMatrix_num(test, rMode, tStep);

  
}



void
NlBeam_SM :: giveInternalForcesVector_u(FloatArray &answer, TimeStep *tStep, const FloatArray &u)
{
  // solution vector
  FloatArray fint(3);
  answer.resize(6);
  FloatArray f_a(3);
  // only the first three entries of f are computed
  this->findLeftEndForces(u, fint);
  answer.at(1) = fint.at(1);
  answer.at(2) = fint.at(2);
  answer.at(3) = fint.at(3);  
  answer.at(4) = -answer.at(1);
  answer.at(5) = -answer.at(2);
  double c1 = beamLength*sin(pitch) + u.at(5) - u.at(2);
  double c2 = -beamLength*cos(pitch) - u.at(4) + u.at(1);
  answer.at(6) = c1 * answer.at(1) + c2 * answer.at(2) - answer.at(3);
}


void
NlBeam_SM :: computeStiffnessMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
  
  answer.resize(6,6);
  FloatArray fint(3), fint_p(3);
  // solution vector
  FloatArray u;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, u);
  this->giveInternalForcesVector_u(fint, tStep, u);   
  double pert = 1.e-6;
  FloatArray up;
  for(int i = 1; i <= 6; i++) {
    up = u;
    up.at(i) += pert;
    this->giveInternalForcesVector_u(fint_p, tStep, up);
    for(int j = 1; j <= 6; j++) {
      answer.at(i,j)  = (fint_p.at(j) - fint.at(j)) / pert;
    }
      
  }

}




Interface*
NlBeam_SM :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case VTKXMLExportModuleElementInterfaceType:
      return static_cast< VTKXMLExportModuleElementInterface * >( this );
    default:
      return StructuralElement :: giveInterface(it);
    }
}



void
NlBeam_SM :: giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
    vtkPieces.resize(1);
 
    int numCells = this->NIP;
    const int numCellNodes  = 2; // linear line
    int nNodes = numCellNodes * numCells;

    vtkPieces.at(0).setNumberOfCells(numCells);
    vtkPieces.at(0).setNumberOfNodes(nNodes);

    int val    = 1;
    int offset = 0;
    IntArray nodes(numCellNodes);
 
    Node *nodeA = this->giveNode(1);
    
    FloatArray uab, ug(3);
    this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
    int nodeNum = 1;
    double ds = beamLength/NIP;
    FloatMatrix T;
    this->construct_T(T, uab.at(3));
    FloatArray nodeCoords(3);
    IntArray connectivity(2);
    for ( int iElement = 1; iElement <= numCells; iElement++ ) {
      for (int iNode = 1; iNode <= numCellNodes; iNode++) {
	double L = (iElement - 1) * ds + (iNode - 1) * ds;
	nodeCoords.at(1) = nodeA->giveCoordinate(1) + L * cos(pitch);
	nodeCoords.at(3) = nodeA->giveCoordinate(3) + L * sin(pitch);
	nodeCoords.at(2) = 0;
	vtkPieces.at(0).setNodeCoords(nodeNum, nodeCoords);
	nodeNum++;
	connectivity.at(iNode) = val++;
      }
      vtkPieces.at(0).setConnectivity(iElement, connectivity);
      offset += 2;
      vtkPieces.at(0).setOffset(iElement, offset);
      vtkPieces.at(0).setCellType(iElement, 3);
    }


    int n = primaryVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfPrimaryVarsToExport(n, nNodes);
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = ( UnknownType ) primaryVarsToExport.at(i);
        if ( utype == DisplacementVector ) {
	  FloatArray l;
	  FloatMatrix T;
	  for ( int nN = 1; nN <= nNodes; nN++ ) {
	    int lN = nN % 2;
	    int iNode;
	    if(lN == 0) {
	      iNode = nN/2 + 1;
	    } else {
	      iNode = (nN + 1) / 2;
	    }

	    double L = (iNode-1) * ds;
	    this->construct_l(l, uab.at(3), L);
	    this->construct_T(T, uab.at(3));
	    FloatArray u_l, u_g;
	    u_l = {this->u.at(iNode), this->w.at(iNode), 0};
	    u_l.subtract(l);
	    u_g.beTProductOf(T, u_l);
	    ug.at(1) = u_g.at(1) + uab.at(1);
	    ug.at(3) = u_g.at(2) + uab.at(2);
	    ug.at(2) = 0;
	    //rotation is not exported this->phi.at(i) + uab.at(3);
	    vtkPieces.at(0).setPrimaryVarInNode(i, nN, ug);
	  }
        }
    }


    InternalStateType isttype;
    n = internalVarsToExport.giveSize();
    vtkPieces [ 0 ].setNumberOfInternalVarsToExport(n, nNodes);
    for ( int i = 1; i <= n; i++ ) {
        isttype = ( InternalStateType ) internalVarsToExport.at(i);
	if ( isttype == IST_BeamForceMomentumTensor ) {
	   FloatArray endForces(3);
	   for ( int nN = 1; nN <= nNodes; nN++ ) {
	     int lN = nN % 2;
	     int iNode;
	     if(lN == 0) {
	       iNode = nN/2 + 1;
	     } else {
	       iNode = (nN + 1) / 2;
	     }
	     endForces.at(1) = this->vN.at(iNode);
	     endForces.at(3) = this->vM.at(iNode);
	     vtkPieces [ 0 ].setInternalVarInNode(i, nN, endForces);
	   }
	} else {
	  fprintf( stderr, "VTKXMLExportModule::exportIntVars: unsupported variable type %s\n", __InternalStateTypeToString(isttype) );
	}
    }


    

}



void
NlBeam_SM :: printOutputAt(FILE *file, TimeStep *tStep)
{
  /*
  FILE *FID;
  

  char fext[100];
  sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );
  std :: string fileName, functionname, temp;
  fileName = this->giveDomain()->giveEngngModel()->giveOutputBaseFileName();
  size_t foundDot;
  foundDot = fileName.rfind(".");
  
  while (foundDot != std :: string :: npos) {
    fileName.replace(foundDot, 1, "_");
    foundDot = fileName.rfind(".");   
  }
  
  fileName += fext;
  temp = fileName;
  size_t backslash = temp.rfind("/");
  if (backslash != std :: string :: npos ) {
    functionname = temp.substr(backslash+1, std :: string :: npos);
  } else {
    functionname = temp;
  }


  
  fileName += ".m";
  if ( ( FID = fopen(fileName.c_str(), "w") ) == NULL ) {
    OOFEM_ERROR("failed to open file %s", fileName.c_str() );
  }



  //transform displacements to global coordinate system
  FloatArray uab, ug(NIP+1), wg(NIP+1), phig(NIP+1), u_l, u_g;
  this->computeVectorOf({D_u, D_w, R_v}, VM_Total, tStep, uab);
  ug.at(1) = uab.at(1);
  wg.at(1) = uab.at(2);
  phig.at(1) = uab.at(3);
  double L = 0;
  double dx = beamLength/NIP;
  for (int i=2; i <= NIP+1; i++) {
    L = L + dx;
    FloatArray l;
    FloatMatrix T;
    this->construct_T(T, uab.at(3));
    this->construct_l(l, uab.at(3), L);
    u_l = {this->u.at(i), this->w.at(i), 0};
    u_l.subtract(l);
    u_g.beTProductOf(T, u_l);
    ug.at(i) = u_g.at(1) + ug.at(1);
    wg.at(i) = u_g.at(2) + wg.at(1);
    phig.at(i) = this->phi.at(i) + phig.at(1);
    
  }

  
  

  fprintf( FID, " function [x z u w phi N M]= %s \n", functionname.c_str() );

   fprintf(FID, "x=[");
   for ( double val: x ) {
     fprintf( FID, "%f,", val * cos(pitch) );
   }   
   fprintf(FID, "];\n");
   
   fprintf(FID, "z=[");
   for ( double val: x ) {
     fprintf( FID, "%f,", val * sin(pitch) );
   }   
   fprintf(FID, "];\n");

   
   
   fprintf(FID, "u=[");
   for ( double val: ug ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");
   
   
   fprintf(FID, "w=[");
   for ( double val: wg ) {
     fprintf( FID, "%f,", val );
   }     
   fprintf(FID, "];\n");
   
   fprintf(FID, "phi=[");
   for ( double val: phig ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");


   fprintf(FID, "N=[");
   for ( double val: vN ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");


   fprintf(FID, "M=[");
   for ( double val: vM ) {
     fprintf( FID, "%f,", val );
   }
   fprintf(FID, "];\n");

   
   fprintf(FID, "end\n");

   fclose(FID);
  */
}
  

FILE *
NlBeam_SM :: giveOutputStream(TimeStep *tStep)
{
    FILE *answer;

    char fext[100];
    sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );
    std :: string fileName;
    fileName = this->giveDomain()->giveEngngModel()->giveOutputBaseFileName();
    fileName += fext;
    fileName += ".m";
    if ( ( answer = fopen(fileName.c_str(), "w") ) == NULL ) {
        OOFEM_ERROR("failed to open file %s", fileName.c_str() );
    }
    return answer;
}
  

} // end namespace oofem
