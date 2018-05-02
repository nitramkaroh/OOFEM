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

#include "../sm/EngineeringModels/POD/reducedstate.h"
#include "floatarray.h"
#include "error.h"
#include "datastream.h"

namespace oofem {



IRResultType
ReducedState :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    epsilonPOD = 1.e-9;
    IR_GIVE_OPTIONAL_FIELD(ir, this->epsilonPOD, _IFT_ReducedState_tolortho);

    numberOfStressComponents = 0;
    reducedBasisMatrix.resize(0,0);
    reducedBasisMatrix.zero();


    return IRRT_OK;
}
  
  /*
void 
ReducedState :: subspaceExpansion(FloatArray &FEM_vars, FloatMatrix &reducedBasisMatrix, FloatMatrix &rbCoords, FloatMatrix &covarianceMatrix, double weight) 
{
  */

void 
ReducedState :: subspaceExpansion(FloatArray &FEM_vars, double weight) 
{
  /*
  FloatArray FEM_vars;
  this->takeSnapshot_dofs(FEM_vars);
  */
   FloatArray residuum;
   residuum = FEM_vars;
   int nVectors = reducedBasisMatrix.giveNumberOfColumns();
   //orthogonalize FEM_var against columns of reducedBasisMatrix   
   this->orthogonalize(residuum, reducedBasisMatrix);   
   // if projection error, normalize projection residual and add it to base
   //@todo check the norms and define eps
   
   if ( fabs(residuum.at(residuum.giveIndexMaxAbsElem())) > 1.e-6*fabs(FEM_vars.at(FEM_vars.giveIndexMaxAbsElem())) ) {
     residuum.times(1./norm(residuum));
     reducedBasisMatrix.resizeWithData(FEM_vars.giveSize(),nVectors+1);
     reducedBasisMatrix.setColumn(residuum,nVectors+1);
   }

   // add related reduced variables to rbCoords, a = A'*FEM_var
   FloatArray a;
   a.beTProductOf(reducedBasisMatrix, FEM_vars);
   // resize reduced coordinates matrix
   rbCoords.resizeWithData(a.giveSize(),rbCoords.giveNumberOfColumns()+1);
   rbCoords.setColumn(a,rbCoords.giveNumberOfColumns());
   
   // update covariance:  Covariance += weight*(a*a')
   if(a.giveSize()>covarianceMatrix.giveNumberOfRows()) {
     covarianceMatrix.resizeWithData(a.giveSize(),a.giveSize());
   }
   FloatMatrix aa;
   aa.beDyadicProductOf(a,a);
   covarianceMatrix.add(aa);
   
}


void 
ReducedState :: subspaceExpansion_dofWeights(FloatArray &FEM_vars, FloatArray &dofWeights, double weight) 
{
   FloatArray residuum;
   residuum = FEM_vars;
   int nVectors = reducedBasisMatrix.giveNumberOfColumns();
   //orthogonalize FEM_var against columns of reducedBasisMatrix   
   //   this->orthogonalize(residuum, reducedBasisMatrix);
   this->orthogonalize(residuum, reducedBasisMatrix, dofWeights);

   FloatArray wFEM_vars, wResiduum;
   wFEM_vars = FEM_vars;
   wResiduum = residuum;
   for(int i = 1; i <= dofWeights.giveSize(); i ++) {
     wFEM_vars.at(i) *= dofWeights.at(i);
     wResiduum.at(i) *= dofWeights.at(i);   
   }

   // if projection error, normalize projection residual and add it to base
   //@todo check the norms and define eps
   
   if ( fabs(wResiduum.at(wResiduum.giveIndexMaxAbsElem())) > 1.e-6*fabs(wFEM_vars.at(wFEM_vars.giveIndexMaxAbsElem())) ) {
     for (int i = 1; i <= dofWeights.giveSize(); i++) {
       wResiduum.at(i) *= dofWeights.at(i);
     }
     wResiduum.times(1./norm(wResiduum));
     reducedBasisMatrix.resizeWithData(FEM_vars.giveSize(),nVectors+1);
     reducedBasisMatrix.setColumn(wResiduum,nVectors+1);
   }


   //@todo checkt this
    for(int i = 1; i <= dofWeights.giveSize(); i ++) {
     wFEM_vars.at(i) *= dofWeights.at(i);
     }
   // add related reduced variables to rbCoords, a = A'*FEM_var
   FloatArray a;
   a.beTProductOf(reducedBasisMatrix, wFEM_vars);
   // resize reduced coordinates matrix
   rbCoords.resizeWithData(a.giveSize(),rbCoords.giveNumberOfColumns()+1);
   rbCoords.setColumn(a,rbCoords.giveNumberOfColumns());
   
   // update covariance:  Covariance += weight*(a*a')
   if(a.giveSize()>covarianceMatrix.giveNumberOfRows()) {
     covarianceMatrix.resizeWithData(a.giveSize(),a.giveSize());
   }
   FloatMatrix aa;
   aa.beDyadicProductOf(a,a);
   covarianceMatrix.add(aa);
}


  /*
void 
ReducedState :: subspaceExpansion_dofs(FloatArray &FEM_vars, double weight)
{
  this->subspaceExpansion(FEM_vars, reducedBasisMatrix_dofs, rbCoords_dofs, covarianceMatrix_dofs, weight);
}

void 
ReducedState :: subspaceExpansion_stress(FloatArray &FEM_vars, double weight)
{
  this->subspaceExpansion(FEM_vars, reducedBasisMatrix_stress, rbCoords_stress, covarianceMatrix_stress, weight);
}
  */


  /*bool 
ReducedState :: subspaceSelection(FloatMatrix &reducedBasisMatrix, FloatMatrix &rbCoords, FloatMatrix &covarianceMatrix)
{*/
bool
ReducedState :: subspaceSelection()
{

   FloatArray W;
   FloatMatrix U, Ucut;
   // compute [U,W] = eigs(U); 
   // on exit U hold eigenvectors, W eigenvalues sorted in descending order
   int info = covarianceMatrix.jaco_(W, U, 12);
   if (info !=0) {
     OOFEM_ERROR("Error in computation of eigenvaluse");
   }
   // find the most relevant eigenvalues
   double treshold = epsilonPOD * W.at(W.giveIndexMaxElem());
   int cutoff;
   for (cutoff = 0; cutoff < W.giveSize(); cutoff++) {
     if (W(cutoff) < treshold) { 
	break;
     }
   }
   // restriction of eigenvectors to selected ones, U = U(:,1:cutoff);
   Ucut.beSubMatrixOf(U,1,U.giveNumberOfRows(), 1, cutoff);
    FloatMatrix tmpMatrix;
   // reducedBasisMatrix update : reducedBasisMatrix = reducedBasisMatrix*U;
   tmpMatrix.beProductOf(reducedBasisMatrix,Ucut);
   reducedBasisMatrix = tmpMatrix;
   // Covariance update : Covariance =  U'*Covariance*U;
   tmpMatrix.beProductOf(covarianceMatrix, Ucut);
   covarianceMatrix.beTProductOf(Ucut, tmpMatrix);
   // rbCoords update : rbCoords = U'*rbCoords;
    tmpMatrix.beTProductOf(U, rbCoords);
    rbCoords = tmpMatrix;
  
   return true;
}

  /*
bool 
ReducedState :: subspaceSelection_dofs()
{
  return this->subspaceSelection(reducedBasisMatrix_dofs, rbCoords_dofs, covarianceMatrix_dofs);
}


bool 
ReducedState :: subspaceSelection_stress()
{
  return this->subspaceSelection(reducedBasisMatrix_stress, rbCoords_stress, covarianceMatrix_stress);
}
  */
  

void 
ReducedState :: orthogonalize(FloatArray &answer, const FloatMatrix& A)
{

  for ( int k=1; k<=A.giveNumberOfColumns(); k++) {
    // compute d{k} = -1.*A(:,k)'*R{k}
    double d;
    FloatArray v;
    v.beColumnOf(A,k);
    d = - dot(v, answer);
    // compute R{k+1} = R{k}+d{k}*A(:,k)
    answer.add(d,v);
  }
}


void 
ReducedState :: orthogonalize(FloatArray &answer, const FloatMatrix& A, const FloatArray &dofWeightsArray)
{


  FloatArray Mv = answer;
  for ( int k=1; k<=A.giveNumberOfColumns(); k++) {
    for(int i = 1; i <= dofWeightsArray.giveSize(); i++) {
      Mv.at(i) *= dofWeightsArray.at(i) * dofWeightsArray.at(i); 
    }
    // compute d{k} = -1.*A(:,k)'*R{k}
    double d;
    FloatArray v;
    v.beColumnOf(A,k);
    d = - dot(v, answer);
    // compute R{k+1} = R{k}+d{k}*A(:,k)
    answer.add(d,v);
  }
}



void
ReducedState :: storeYourself(DataStream &stream) 
{
  rbCoords.storeYourself(stream);
  reducedBasisMatrix.storeYourself(stream);
  covarianceMatrix.storeYourself(stream);
  stream.write(numberOfStressComponents);
}


void
ReducedState :: restoreYourself(DataStream &stream) 
{
  rbCoords.restoreYourself(stream);
  reducedBasisMatrix.restoreYourself(stream);
  covarianceMatrix.restoreYourself(stream);
  stream.read(numberOfStressComponents);
}



} // end namespace oofem
