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

#include "lsmastermat.h"
#include "Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(LargeStrainMasterMaterial);

// constructor
LargeStrainMasterMaterial :: LargeStrainMasterMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{
    slaveMat = 0;
}

// destructor
LargeStrainMasterMaterial :: ~LargeStrainMasterMaterial()
{ }

// reads the model parameters from the input file
IRResultType
LargeStrainMasterMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    IR_GIVE_OPTIONAL_FIELD(ir, slaveMat, _IFT_LargeStrainMasterMaterial_slaveMat); // number of slave material
    IR_GIVE_OPTIONAL_FIELD(ir, m, _IFT_LargeStrainMasterMaterial_m); // type of Set-Hill strain tensor

    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
LargeStrainMasterMaterial :: CreateStatus(GaussPoint *gp) const
{
  return new LargeStrainMasterMaterialStatus(1, this->giveDomain(), gp, slaveMat);
}


void
LargeStrainMasterMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
{
    LargeStrainMasterMaterialStatus *status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);

    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );

    FloatArray eVals, vSethHillStrain, vSethHillStress, lambda, eps(3);
    FloatMatrix F, C, eVecs, SethHillStrain, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(vF);
    //compute right Cauchy-Green tensor(C)
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    // compute Seth - Hill's strain measure, it depends on parameter m
    lambda = eVals;
    if ( m == 0 ) {
      eps.at(1) = 1. / 2. * log(lambda.at(1));
      eps.at(2) = 1. / 2. * log(lambda.at(2));
      eps.at(3) = 1. / 2. * log(lambda.at(3));
    } else {
      eps.at(1) = 1. /  m  * ( pow(lambda.at(1), m / 2. ) - 1. );
      eps.at(2) = 1. /  m  * ( pow(lambda.at(2), m / 2. ) - 1. );
      eps.at(3) = 1. /  m  * ( pow(lambda.at(3), m / 2. ) - 1. );
    }

    SethHillStrain.resize(3, 3);
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
	  SethHillStrain.at(i, j) = eps.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + eps.at(2) *eVecs.at(i, 2) * eVecs.at(j, 2) + eps.at(3) *eVecs.at(i, 3) * eVecs.at(j, 3);  
        }
    }

    vSethHillStrain.beSymVectorFormOfStrain(SethHillStrain);
    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->giveRealStressVector_3d(vSethHillStress, slaveGp, vSethHillStrain, tStep);
    
    SethHillStress.beMatrixForm(vSethHillStress);
    // transformation matrices
    this->giveTransformationMatrices(PP,TL, F, SethHillStress, eVals, eVecs);
    answer.beTProductOf(PP, vSethHillStress);

    status->letTempStressVectorBe(vSethHillStress);
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);

}


void
LargeStrainMasterMaterial :: giveCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
{
    LargeStrainMasterMaterialStatus *status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );
    this->initTempStatus(gp);

    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );

    FloatArray eVals, vSethHillStrain, vSethHillStress, lambda, eps(3);
    FloatMatrix F, b, eVecs, SethHillStrain, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(vF);
    //compute right Cauchy-Green tensor(C)
    b.beProductTOf(F, F);
    // compute eigen values and eigen vectors of C
    b.jaco_(eVals, eVecs, 15);
    // compute Seth - Hill's strain measure, it depends on parameter m
    lambda = eVals;
    if ( m == 0 ) {
      eps.at(1) = 1. / 2. * log(lambda.at(1));
      eps.at(2) = 1. / 2. * log(lambda.at(2));
      eps.at(3) = 1. / 2. * log(lambda.at(3));
    } else {
      eps.at(1) = 1. /  m  * ( pow(lambda.at(1), m / 2. ) - 1. );
      eps.at(2) = 1. /  m  * ( pow(lambda.at(2), m / 2. ) - 1. );
      eps.at(3) = 1. /  m  * ( pow(lambda.at(3), m / 2. ) - 1. );
    }

    SethHillStrain.resize(3, 3);
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
	  SethHillStrain.at(i, j) = eps.at(1) * eVecs.at(i, 1) * eVecs.at(j, 1) + eps.at(2) *eVecs.at(i, 2) * eVecs.at(j, 2) + eps.at(3) *eVecs.at(i, 3) * eVecs.at(j, 3);
        }
    }

    vSethHillStrain.beSymVectorFormOfStrain(SethHillStrain);
    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->giveRealStressVector_3d(vSethHillStress, slaveGp, vSethHillStrain, tStep);
    
    SethHillStress.beMatrixForm(vSethHillStress);
    // transformation matrices
    this->giveTransformationMatrices(PP,TL, F, SethHillStress, eVals, eVecs);
    answer.beTProductOf(PP, vSethHillStress);
    // store F, P, and Seth-Hill stress into status
    status->letTempStressVectorBe(vSethHillStress);
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
}



void
LargeStrainMasterMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    LargeStrainMasterMaterialStatus *status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    FloatArray eVals;
    FloatMatrix stiffness, F, C, P, eVecs, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(status->giveTempFVector());
    //compute right Cauchy-Green tensor(C), its eigenvalues and eigenvectors
    C.beTProductOf(F, F);
    // compute eigen values and eigen vectors of C
    C.jaco_(eVals, eVecs, 15);
    

    FloatArray vSethHillStress, vP(status->giveTempPVector());

    if(gp->giveMaterialMode() != _3dMat) {
      StructuralMaterial :: giveFullSymVectorForm( vSethHillStress, status->giveTempStressVector(), gp->giveMaterialMode() );
    } else {
      vSethHillStress =  status->giveTempStressVector();
    }


    P.beMatrixForm(vP);    
    SethHillStress.beMatrixForm(vSethHillStress);
    FloatMatrix invF, S, EP, delta_S;;
    // compute 2-PK stress
    invF.beInverseOf(F);
    S.beProductOf(invF, P);

    this->giveTransformationMatrices(PP,TL, F, SethHillStress, eVals, eVecs);
    this->giveDeltaS_Product(delta_S, S);


    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->give3dMaterialStiffnessMatrix(stiffness, mode, slaveGp, tStep);
    EP.beProductOf(stiffness, PP);

    answer.beTProductOf(PP, EP);
    answer.add(TL);   
    answer.add(delta_S);
    

}


void
LargeStrainMasterMaterial :: giveSpatial3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    LargeStrainMasterMaterialStatus *status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    FloatArray eVals;
    FloatMatrix stiffness, F, b, C, eVecs, SethHillStress, PP, TL;

    //store of deformation gradient into 3x3 matrix
    F.beMatrixForm(status->giveTempFVector());
    //compute right Cauchy-Green tensor(C), its eigenvalues and eigenvectors
    b.beProductTOf(F, F);
    // compute eigen values and eigen vectors of C
    b.jaco_(eVals, eVecs, 15);
    

    FloatArray vSethHillStress, vC(status->giveTempCVector());

    if(gp->giveMaterialMode() != _3dMat) {
      StructuralMaterial :: giveFullSymVectorForm( vSethHillStress, status->giveTempStressVector(), gp->giveMaterialMode() );
    } else {
      vSethHillStress =  status->giveTempStressVector();
    }

    C.beMatrixForm(vC);    
    SethHillStress.beMatrixForm(vSethHillStress);
    FloatMatrix invF, S, EP, delta_S;;

    this->giveTransformationMatrices(PP,TL, F, S, eVals, eVecs);
    this->giveDeltaS_Product(delta_S, S);


    GaussPoint *slaveGp = status->giveSlaveGp();
    sMat->give3dMaterialStiffnessMatrix(stiffness, mode, slaveGp, tStep);
    EP.beProductOf(stiffness, PP);

    answer.beTProductOf(PP, EP);
    answer.add(TL);   
    answer.add(delta_S);
    
}





void
LargeStrainMasterMaterial :: giveDeltaS_Product(FloatMatrix &answer, const FloatMatrix &S)
{   

    answer.resize(9,9);
  
    answer.at(1,1) = S.at(1,1);
    answer.at(1,5) = S.at(1,3);
    answer.at(1,6) = S.at(1,2);
    
    answer.at(2,2) = S.at(2,2);
    answer.at(2,4) = S.at(2,3);
    answer.at(2,9) = S.at(1,2);
    
    answer.at(3,3) = S.at(3,3);
    answer.at(3,7) = S.at(2,3);
    answer.at(3,8) = S.at(2,3);
    
    
    answer.at(4, 2) = S.at(3,2);
    answer.at(4, 4) = S.at(3,3);
    answer.at(4, 9) = S.at(3,1);
    
    answer.at(5, 1) = S.at(3,1);
    answer.at(5, 5) = S.at(3,3);
    answer.at(5, 6) = S.at(3,2);
    
    answer.at(6, 1) = S.at(2,1);
    answer.at(6, 5) = S.at(2,3);
    answer.at(6, 6) = S.at(2,2);
    
    answer.at(7, 3) = S.at(2,3);
    answer.at(7, 7) = S.at(2,2);
    answer.at(7, 8) = S.at(2,1);
    
    answer.at(8, 3) = S.at(1,3);
    answer.at(8, 7) = S.at(1,2);
    answer.at(8, 8) = S.at(1,1);
    
    answer.at(9, 2) = S.at(1,2);
    answer.at(9, 4) = S.at(1,3);
    answer.at(9, 9) = S.at(1,1);
}



void
LargeStrainMasterMaterial :: giveTransformationMatrices(FloatMatrix &PP,FloatMatrix &TL, const FloatMatrix &F, const FloatMatrix &SetHillStress, const FloatArray &lam, const FloatMatrix &N)
{   

    FloatMatrix n;
    //n = N;
    n.beProductOf(F,N);
    FloatArray  d(3),f(3), eps(3);
    if(this->m == 0) { // log strain formulation (m = 0) has a special treatment
      d.at(1)= 1./lam.at(1);
      d.at(2)= 1./lam.at(2);
      d.at(3)= 1./lam.at(3);
      
      f.at(1) = -2./lam.at(1)/lam.at(1);
      f.at(2) = -2./lam.at(2)/lam.at(2);
      f.at(3) = -2./lam.at(3)/lam.at(3);

      eps.at(1) = log(lam.at(1))/2.;
      eps.at(2) = log(lam.at(2))/2.;
      eps.at(3) = log(lam.at(3))/2.;
    } else {  // all other strain measures (m != 0)
      d.at(1)= pow(lam.at(1), m/2.-1);
      d.at(2)= pow(lam.at(2), m/2.-1);
      d.at(3)= pow(lam.at(3), m/2.-1);
      
      f.at(1) = 2. * ( m/2. - 1 ) * pow( lam.at(1), m/2. - 2. );
      f.at(2) = 2. * ( m/2. - 1 ) * pow( lam.at(2), m/2. - 2. );
      f.at(3) = 2. * ( m/2. - 1 ) * pow( lam.at(3), m/2. - 2. );
      
      eps.at(1) = 1./m * ( pow(lam.at(1),m/2.) - 1. );
      eps.at(2) = 1./m * ( pow(lam.at(2),m/2.) - 1. );
      eps.at(3) = 1./m * ( pow(lam.at(3),m/2.) - 1. );

    }

    double eta = 0.;
    FloatMatrix TN, dzeta(3,3), ksi(3,3), theta(3,3);

    TN.beProductOf(SetHillStress,N);
    dzeta.beTProductOf(N, TN);
    




    // compute auxiliary variables 
    // the computation differes depends on if the eigenvalues of C are equal or not
    if(lam.at(1) != lam.at(2)) {
      if(lam.at(2) != lam.at(3)) {
	if(lam.at(1) != lam.at(3)) {
	  // all eigenvalues are different
	  for(int i = 1; i <= 3; i++) {
	    for(int j = 1; j <= 3; j++) {
	      if(i == j) {
		continue;
	      } else {
		theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
		ksi.at(i,j) = (theta.at(i,j) - 1./2.*d.at(j))/(lam.at(i) - lam.at(j));
		for(int k = 1; k <= 3; k++) {
		  if((k != i) && (k != j)) {
		    eta += eps.at(i)/2./(lam.at(i) - lam.at(j))/(lam.at(i)-lam.at(k));
		  }
		}
	      }
	    }
	  }
	} else { //l1 == l3 && l1 != l2
	  for(int i = 1; i <= 3; i++) {
	    for(int j = 1; j <= 3; j++) {
	      if( i == j ) {
		continue;
	      } else {
		if((i == 1 && j == 3) || (i == 3 && j == 1) ) {
		  theta.at(i,j) = 1./2.*d.at(i);
		  ksi.at(i,j) = 1./8.*f.at(i);
		} else {
		  theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
		  ksi.at(i,j) = (theta.at(i,j) - 1./2.*d.at(j))/(lam.at(i) - lam.at(j));
		}
	      }
	    }
	  }	  
	  eta = ksi.at(1,2);
	}
      } else { //l2 == l3 && l1 != l2
	for(int i = 1; i <= 3; i++) {
	  for(int j = 1; j <= 3; j++) {
	    if( i == j ) {
	      continue;
	    } else {
	      if((i == 2 && j == 3) || (i == 3 && j == 2) ) {
		theta.at(i,j) = 1./2.*d.at(i);
		ksi.at(i,j) = 1./8.*f.at(i);
	      } else {
		theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
		ksi.at(i,j) = (theta.at(i,j) - 1./2.*d.at(j))/(lam.at(i) - lam.at(j));
	      }
	    }
	  }
	} 
	eta = ksi.at(1,2);    
      }
    } else if(lam.at(1) != lam.at(3)) { // l1 == l2  && l1 != l3
      for(int i = 1; i <= 3; i++) {
	for(int j = 1; j <= 3; j++) {
	  if( i == j ) {
	    continue;
	  } else {
	    if((i == 1 && j == 2) || (i == 2 && j == 1) ) {
	      theta.at(i,j) = 1./2.*d.at(i);
	      ksi.at(i,j) = 1./8.*f.at(i);
	    } else {
	      theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
	      ksi.at(i,j) = (theta.at(i,j) - 1./.2*d.at(j))/(lam.at(i) - lam.at(j));
	    }
	  }
	}
      } 
      eta = ksi.at(1,3);
    } else {  // l1 == l2 == l3
      for(int i = 1; i <= 3; i++) {
	for(int j = 1; j <= 3; j++) {
	  theta.at(i,j) = 1./2. * d.at(i);
	  ksi.at(i,j) = 1./8. * f.at(i);
	}
      }
      eta = 1./8. * f.at(1);
    }


  FloatMatrix M(9,9);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
	for (int l = 1; l <=3; l++) {
	  M.at(giveVI(i,j),giveVI(k,l)) = n.at(k,i)*N.at(l,j) + n.at(k,j)*N.at(l,i);
	}
      }
    }
  }
    
  PP.resize(6,9), TL.resize(9,9);
  for (int k = 1; k <= 3; k++) {
    for (int l = 1; l <= 3; l++) {
      for (int m = 1; m <=3; m++) {
	for (int n = 1; n<=3; n++) {
	  for (int i = 1; i <= 3; i++) {
	    if(k == l ) {
	      PP.at(giveVI(k,l),giveVI(m,n)) += 0.5 * d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    } else if(k < l) {
	      PP.at(giveVI(k,l),giveVI(m,n)) += d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    }

	    TL.at(giveVI(k,l),giveVI(m,n)) += 1./4.*f.at(i)*dzeta.at(i,i)*M.at(giveVI(i,i),giveVI(k,l))*M.at(giveVI(i,i),giveVI(m,n));
	    for (int j  = 1; j <= 3; j++) {
	      if(j != i) {
		if( k == l) {
		  PP.at(giveVI(k,l),giveVI(m,n)) += theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
		} else if(k < l) {
		  PP.at(giveVI(k,l),giveVI(m,n)) += 2. * theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
		}
		TL.at(giveVI(k,l),giveVI(m,n)) += 2.*ksi.at(i,j) * ( dzeta.at(i,j)*(M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(j,j),giveVI(m,n)) + M.at(giveVI(j,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n))) + dzeta.at(j,j)*M.at(giveVI(i,j),giveVI(k,l))*M.at(giveVI(i,j),giveVI(m,n) ));
		for (int q = 1; q <= 3; q++) {
		  if(q == i) {
		    continue;
		  } else if(q == j) {
		    continue;
		  } else {
		    TL.at(giveVI(k,l),giveVI(m,n)) += 2.*eta*dzeta.at(i,j)*M.at(giveVI(i,q),giveVI(k,l))*M.at(giveVI(j,q),giveVI(m,n));
		  }
		}
	      }
	    }
	  }
	}	    
      }
    }
  }   



}


int
LargeStrainMasterMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    LargeStrainMasterMaterialStatus *status = static_cast< LargeStrainMasterMaterialStatus * >( this->giveStatus(gp) );


    GaussPoint *slaveGp = status->giveSlaveGp();
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    return sMat->giveIPValue(answer, slaveGp, type, tStep);
}


//=============================================================================

  LargeStrainMasterMaterialStatus :: LargeStrainMasterMaterialStatus(int n, Domain *d, GaussPoint *g, int s) : StructuralMaterialStatus(n, d, g),
    Pmatrix(6, 6),
    TLmatrix(6, 6),
    transformationMatrix(6, 6),
    slaveMat(s)
{
    Pmatrix.beUnitMatrix();
    FloatArray coords = {0,0,0};
    slaveGp = new GaussPoint(NULL, 1, coords,1.0, _3dMat );
    
}


LargeStrainMasterMaterialStatus :: ~LargeStrainMasterMaterialStatus()
{ }


void
LargeStrainMasterMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    MaterialStatus *mS = sMat->giveStatus(slaveGp);

    mS->printOutputAt(file, tStep);
    //  StructuralMaterialStatus :: printOutputAt(file, tStep);
}


// initializes temporary variables based on their values at the previous equlibrium state
void LargeStrainMasterMaterialStatus :: initTempStatus()
{
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    MaterialStatus *mS = sMat->giveStatus(slaveGp);
    mS->initTempStatus();

    if ( this->giveCVector().giveSize() == 0 ) {
      stressVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
    }

    StructuralMaterialStatus :: initTempStatus();
    tempCVector      = CVector;
      
    //
}


// updates internal variables when equilibrium is reached
void
LargeStrainMasterMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterial *sMat = static_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    MaterialStatus *mS = sMat->giveStatus(slaveGp);
    mS->updateYourself(tStep);
    StructuralMaterialStatus :: updateYourself(tStep);
    CVector      = tempCVector;
	
    
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
LargeStrainMasterMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    StructuralMaterial *sMat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial(slaveMat) );
    MaterialStatus *mS = sMat->giveStatus(slaveGp);
    // save parent class status
    if ( ( iores = mS->saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    return CIO_OK;
}


contextIOResultType
LargeStrainMasterMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem




