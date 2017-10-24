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

#include "misesmatlogstrain.h"
#include "Materials/misesmat.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(MisesMatLogStrain);

// constructor
  MisesMatLogStrain :: MisesMatLogStrain(int n, Domain *d) : MisesMat(n, d)
{

}

// destructor
MisesMatLogStrain :: ~MisesMatLogStrain()
{
}


// creates a new material status  corresponding to this class
MaterialStatus *
MisesMatLogStrain :: CreateStatus(GaussPoint *gp) const
{
    return new MisesMatLogStrainStatus(1, this->giveDomain(), gp);
}



void
MisesMatLogStrain :: giveCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vf, TimeStep *tStep)
{
   // initialization
    this->initTempStatus(gp);
    MisesMatLogStrainStatus *status = static_cast< MisesMatLogStrainStatus * >( this->giveStatus(gp) );
    double J;
    FloatArray vFn, vF, eps_n, e(3), vEps, lambda, plStrain;
    FloatMatrix Fn, f, F, Bn, Bf, Btr, eps, N, mEps;
    
    //store of deformation gradient into 3x3 matrix
    vFn = status->giveFVector();
    Fn.beMatrixForm(vFn);
    f.beMatrixForm(vf);
    F.beProductOf(f,Fn);
    vF.beVectorForm(F);
    J = F.giveDeterminant();
    //computes elastic trial left Cauchy-Green tensor(B trial)
    eps_n = status->giveStrainVector();
    plStrain = status->givePlasticStrain();
    eps_n.subtract(plStrain);
    mEps.beMatrixFormOfStrain(eps_n);
    mEps.jaco_(lambda, N, 15);
    FloatArray eLambda(3);
    eLambda.at(1) = exp(2.*lambda.at(1));
    eLambda.at(2) = exp(2.*lambda.at(2));
    eLambda.at(3) = exp(2.*lambda.at(3));
    Bn.resize(3,3);
    for ( int i = 1; i < 4; i++ ) {
      for ( int j = 1; j < 4; j++ ) {
	Bn.at(i, j) = eLambda.at(1) * N.at(i, 1) * N.at(j, 1) + eLambda.at(2) * N.at(i, 2) * N.at(j, 2) + eLambda.at(3) * N.at(i, 3) * N.at(j, 3);
      }
    }
    Bf.beProductTOf(Bn,f);
    Btr.beProductOf(f,Bf);
    //    Btr.beProductTOf(F,F);
    // compute eigen values and eigen vectors of B trial
    Btr.jaco_(lambda, N, 15);
    // compute logarithmic strain
    e.at(1) = 1. / 2. * log(lambda.at(1));
    e.at(2) = 1. / 2. * log(lambda.at(2));
    e.at(3) = 1. / 2. * log(lambda.at(3));
    eps.resize(3, 3);
    for ( int i = 1; i < 4; i++ ) {
        for ( int j = 1; j < 4; j++ ) {
	  eps.at(i, j) = e.at(1) * N.at(i, 1) * N.at(j, 1) + e.at(2) * N.at(i, 2) * N.at(j, 2) + e.at(3) * N.at(i, 3) * N.at(j, 3);
	}
    }
    // add plastic strain to obtaine total strain
    vEps.beSymVectorFormOfStrain(eps);
    vEps.add(plStrain);
    // compute Kirchhoff stress using the small strain algorithm
    this->giveRealStressVector_3d(answer, gp, vEps, tStep);
    // transform Kirchhoff stress to Cauchy stress
    answer.times(1./J);
    // store the deformation gradient and Cauchy stress
    status->letTempFVectorBe(vF);
    status->letTempCVectorBe(answer);
}

void
MisesMatLogStrain :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet

{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    double J;
    FloatArray vC, vFn, vf;
    FloatMatrix Fn, F, f, C, invFn, invF, P;
  
    vFn = status->giveFVector();
    Fn.beMatrixForm(vFn);
    invFn.beInverseOf(Fn);
    
    F.beMatrixForm(vF);
    invF.beInverseOf(F);
    
    f.beProductOf(F,invFn);
    vf.beVectorForm(f);

    // compute Cauchy stress
    this->giveCauchyStressVector_3d(vC, gp, vf, tStep);
    C.beMatrixForm(vC);
    // transform Cauchy stress to first PK stress
    J = F.giveDeterminant();  
    P.beProductTOf(C,invF);
    P.times(J);    
    answer.beVectorForm(P);

    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(answer);
}



  

void 
MisesMatLogStrain :: giveSpatial3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{ 
    FloatMatrix stiffness;
    this->give3dMaterialStiffnessMatrix(stiffness, mode, gp, tStep);
    this->convert_LogStiffness_2_SpatialStiffness(answer, stiffness, mode, gp, tStep);
}


void
MisesMatLogStrain :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                           MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *tStep)
// returns the 9x9 tangent stiffness matrix - dP/dF
{

    answer.resize(9,9);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    double J;
    FloatArray vF = status->giveTempFVector();
    FloatMatrix spatialStiffness, F, invF;
    this->giveSpatial3dMaterialStiffnessMatrix(spatialStiffness, mode, gp, tStep);
    
    F.beMatrixForm(vF);
    invF.beInverseOf(F);
    J = F.giveDeterminant();
      
   
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
	for (int k = 1; k <= 3; k++) {
	  for (int l = 1; l <= 3; l++) {
	    for (int m = 1; m <= 3; m++) {
	      for (int n = 1; n <= 3; n++) {
		answer.at(giveVI(i,j),giveVI(k,l)) += J * spatialStiffness.at(giveVI(i,m),giveVI(k,n))*invF.at(j,m)*invF.at(l,n);
	      }
	    }
	  }
	}
      }
    }

   

}




  
void 
MisesMatLogStrain :: convert_LogStiffness_2_SpatialStiffness(FloatMatrix &answer, const FloatMatrix &stiffness,MatResponseMode  mode,  GaussPoint *gp, TimeStep *tStep)
{
    MisesMatLogStrainStatus *status = static_cast< MisesMatLogStrainStatus * >( this->giveStatus(gp) );
    //store of deformation gradient into 3x3 matrix
    double J;
    FloatArray vF, lambda, vCauchy;
    FloatMatrix F, Btr, LL,BB , N, deltaSigma;
    vF = status->giveTempFVector();
    F.beMatrixForm(vF);
    J = F.giveDeterminant();

    // trial elastic strain
    FloatArray vTrialElStrain = status->giveTempStrainVector();
    FloatArray vPlStrain = status->givePlasticStrain();
    vTrialElStrain.subtract(vPlStrain);
    FloatArray vFullTrialElStrain;
    StructuralMaterial :: giveFullSymVectorForm(vFullTrialElStrain, vTrialElStrain, gp->giveMaterialMode());
    FloatMatrix trialElStrain;   
    trialElStrain.beMatrixFormOfStrain(vFullTrialElStrain);
    // spectral decomosition of trial elastic strain
    trialElStrain.jaco_(lambda, N, 15);
    FloatArray eLambda(3);
    eLambda.at(1) = exp(2.*lambda.at(1));
    eLambda.at(2) = exp(2.*lambda.at(2));
    eLambda.at(3) = exp(2.*lambda.at(3));
    Btr.resize(3,3);
    for ( int i = 1; i < 4; i++ ) {
      for ( int j = 1; j < 4; j++ ) {
	Btr.at(i, j) = eLambda.at(1) * N.at(i, 1) * N.at(j, 1) + eLambda.at(2) * N.at(i, 2) * N.at(j, 2) + eLambda.at(3) * N.at(i, 3) * N.at(j, 3);
      }
    }
  

    this->compute_dElog_dE(LL, F, eLambda, N);
    this->compute_dB_dF(BB, Btr);       
    this->giveDLBProduct(answer,stiffness, LL, BB);
    
    answer.times(1./2./J);
    vCauchy = status->giveTempCVector();    
    if(vCauchy.giveSize()) {
      this->give_Sigma_Delta_Product(deltaSigma, vCauchy);
      answer.subtract(deltaSigma);
    }
    




}



void
MisesMatLogStrain :: giveDLBProduct(FloatMatrix &answer, const FloatMatrix &DD, const FloatMatrix &LL, const FloatMatrix &BB)
{

  answer.resize(9,9);


  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
	for (int l = 1; l <= 3; l++) {
	  for (int m = 1; m <= 3; m++) {
	    for (int n = 1; n <= 3; n++) {
	      for (int o = 1; o <= 3; o++) {
		for (int p = 1; p <= 3; p++) {
		  answer.at(giveVI(i,j),giveVI(k,l)) += DD.at(giveSymVI(i,j),giveSymVI(m,n)) * LL.at(giveVI(m,n),giveVI(o,p))*BB.at(giveVI(o,p),giveVI(k,l));
		}
	      }
	    }
	  }
	}
      }
    }
  }
}



void
MisesMatLogStrain :: give_Sigma_Delta_Product( FloatMatrix &answer, const FloatArray &vCauchy)
{
  // product answer_ijkl = sigma_il * delta_jk
  /*[ sig11,     0,     0,     0, sig13, sig12,     0,     0,     0]
    [     0, sig22,     0, sig23,     0,     0,     0,     0, sig21]
    [     0,     0, sig33,     0,     0,     0, sig32, sig31,     0]
    [     0,     0, sig23,     0,     0,     0, sig22, sig21,     0]
    [     0,     0, sig13,     0,     0,     0, sig12, sig11,     0]
    [     0, sig12,     0, sig13,     0,     0,     0,     0, sig11]
    [     0, sig32,     0, sig33,     0,     0,     0,     0, sig31]
    [ sig31,     0,     0,     0, sig33, sig32,     0,     0,     0]
    [ sig21,     0,     0,     0, sig23, sig22,     0,     0,     0]
  */
  answer.resize(9,9);
  answer.at(1,1) = vCauchy.at(1);
  answer.at(1,5) = vCauchy.at(5);
  answer.at(1,6) = vCauchy.at(6);

  answer.at(2,2) = vCauchy.at(2);
  answer.at(2,4) = vCauchy.at(4);
  answer.at(2,9) = vCauchy.at(6);

  answer.at(3,3) = vCauchy.at(3);
  answer.at(3,7) = vCauchy.at(4);
  answer.at(3,8) = vCauchy.at(5);

  answer.at(4,3) = vCauchy.at(4);
  answer.at(4,7) = vCauchy.at(2);
  answer.at(4,8) = vCauchy.at(6);

  answer.at(5,3) = vCauchy.at(5);
  answer.at(5,7) = vCauchy.at(6);
  answer.at(5,8) = vCauchy.at(1);
  
  answer.at(6,2) = vCauchy.at(6);
  answer.at(6,4) = vCauchy.at(5);
  answer.at(6,9) = vCauchy.at(1);

  answer.at(7,2) = vCauchy.at(4);
  answer.at(7,4) = vCauchy.at(3);
  answer.at(7,9) = vCauchy.at(5);

  answer.at(8,1) = vCauchy.at(5);
  answer.at(8,5) = vCauchy.at(3);
  answer.at(8,6) = vCauchy.at(4);

  answer.at(9,1) = vCauchy.at(6);
  answer.at(9,5) = vCauchy.at(4);
  answer.at(9,6) = vCauchy.at(2);

}



void
MisesMatLogStrain :: compute_dB_dF(FloatMatrix &answer, const FloatMatrix &Btr)
{
  // product answer_ijkl = delta_ik * B_jl + delta_jk * B_il
  /*
    [ 2*B11,     0,     0,     0, 2*B13, 2*B12,     0,     0,     0]
    [     0, 2*B22,     0, 2*B23,     0,     0,     0,     0, 2*B21]
    [     0,     0, 2*B33,     0,     0,     0, 2*B32, 2*B31,     0]
    [     0,   B32,   B23,   B33,     0,     0,   B22,   B21,   B31]
    [   B31,     0,   B13,     0,   B33,   B32,   B12,   B11,     0]
    [   B21,   B12,     0,   B13,   B23,   B22,     0,     0,   B11]
    [     0,   B32,   B23,   B33,     0,     0,   B22,   B21,   B31]
    [   B31,     0,   B13,     0,   B33,   B32,   B12,   B11,     0]
    [   B21,   B12,     0,   B13,   B23,   B22,     0,     0,   B11]
  */
  answer.resize(9,9);
  answer.at(1,1) = 2. * Btr.at(1,1);
  answer.at(1,5) = 2. * Btr.at(1,3);
  answer.at(1,6) = 2. * Btr.at(1,2);

  answer.at(2,2) = 2. * Btr.at(2,2);
  answer.at(2,4) = 2. * Btr.at(2,3);
  answer.at(2,9) = 2. * Btr.at(2,1);

  answer.at(3,3) = 2. * Btr.at(3,3);
  answer.at(3,7) = 2. * Btr.at(3,2);
  answer.at(3,8) = 2. * Btr.at(3,1);

  answer.at(4,2) = Btr.at(3,2);
  answer.at(4,3) = Btr.at(2,3);
  answer.at(4,4) = Btr.at(3,3);
  answer.at(4,7) = Btr.at(2,2);
  answer.at(4,8) = Btr.at(2,1);
  answer.at(4,9) = Btr.at(3,1);

  answer.at(5,1) = Btr.at(3,1);
  answer.at(5,3) = Btr.at(1,3);
  answer.at(5,5) = Btr.at(3,3);
  answer.at(5,6) = Btr.at(3,2);
  answer.at(5,7) = Btr.at(1,2);
  answer.at(5,8) = Btr.at(1,1);
  
  answer.at(6,1) = Btr.at(2,1);
  answer.at(6,2) = Btr.at(1,2);
  answer.at(6,4) = Btr.at(1,3);
  answer.at(6,5) = Btr.at(2,3);
  answer.at(6,6) = Btr.at(2,2);
  answer.at(6,9) = Btr.at(1,1);

  answer.at(7,2) = Btr.at(3,2);
  answer.at(7,3) = Btr.at(2,3);
  answer.at(7,4) = Btr.at(3,3);
  answer.at(7,7) = Btr.at(2,2);
  answer.at(7,8) = Btr.at(2,1);
  answer.at(7,9) = Btr.at(3,1);

  answer.at(8,1) = Btr.at(3,1);
  answer.at(8,3) = Btr.at(1,3);
  answer.at(8,5) = Btr.at(3,3);
  answer.at(8,6) = Btr.at(3,2);
  answer.at(8,7) = Btr.at(1,2);
  answer.at(8,8) = Btr.at(1,1);

  answer.at(9,1) = Btr.at(2,1);
  answer.at(9,2) = Btr.at(1,2);
  answer.at(9,4) = Btr.at(1,3);
  answer.at(9,5) = Btr.at(2,3);
  answer.at(9,6) = Btr.at(2,2);
  answer.at(9,9) = Btr.at(1,1);
}


void
MisesMatLogStrain :: compute_dElog_dE(FloatMatrix &answer, const FloatMatrix &F, const FloatArray &lam, const FloatMatrix &N)
{   

    FloatArray  d(3), eps(3);
    FloatMatrix theta(3,3);

    d.at(1)= 1./lam.at(1);
    d.at(2)= 1./lam.at(2);
    d.at(3)= 1./lam.at(3);

    
    eps.at(1) = log(lam.at(1))/2.;
    eps.at(2) = log(lam.at(2))/2.;
    eps.at(3) = log(lam.at(3))/2.;
    
   

    


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
		} else {
		  theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
		}
	      }
	    }
	  }	  
	}
      } else { //l2 == l3 && l1 != l2
	for(int i = 1; i <= 3; i++) {
	  for(int j = 1; j <= 3; j++) {
	    if( i == j ) {
	      continue;
	    } else {
	      if((i == 2 && j == 3) || (i == 3 && j == 2) ) {
		theta.at(i,j) = 1./2.*d.at(i);
	      } else {
		theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
	      }
	    }
	  }
	} 
      }
    } else if(lam.at(1) != lam.at(3)) { // l1 == l2  && l1 != l3
      for(int i = 1; i <= 3; i++) {
	for(int j = 1; j <= 3; j++) {
	  if( i == j ) {
	    continue;
	  } else {
	    if((i == 1 && j == 2) || (i == 2 && j == 1) ) {
	      theta.at(i,j) = 1./2.*d.at(i);
	    } else {
	      theta.at(i,j) = (eps.at(i) - eps.at(j))/(lam.at(i) - lam.at(j));
	    }
	  }
	}
      } 
    } else {  // l1 == l2 == l3
      for(int i = 1; i <= 3; i++) {
	for(int j = 1; j <= 3; j++) {
	  theta.at(i,j) = 1./2. * d.at(i);
	}
      }
    }


  FloatMatrix M(9,9);

  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 3; j++) {
      for (int k = 1; k <= 3; k++) {
	for (int l = 1; l <=3; l++) {
	  M.at(giveVI(i,j),giveVI(k,l)) = N.at(k,i)*N.at(l,j) + N.at(k,j)*N.at(l,i);
	}
      }
    }
  }
    
  answer.resize(9,9);
  for (int k = 1; k <= 3; k++) {
    for (int l = 1; l <= 3; l++) {
      for (int m = 1; m <=3; m++) {
	for (int n = 1; n<=3; n++) {
	  for (int i = 1; i <= 3; i++) {
	      answer.at(giveVI(k,l),giveVI(m,n)) += 0.5 * d.at(i) * N.at(k,i) * N.at(l,i) * M.at(giveVI(i,i),giveVI(m,n));
	    for (int j  = 1; j <= 3; j++) {
	      if(j != i) {
		  answer.at(giveVI(k,l),giveVI(m,n)) += theta.at(i,j)*N.at(k,i)*N.at(l,j)*M.at(giveVI(i,j),giveVI(m,n));
	      }
	    }
	  }
	}	    
      }
    }
  }  
}









int
MisesMatLogStrain :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    MisesMatLogStrainStatus *status = static_cast< MisesMatLogStrainStatus * >( this->giveStatus(gp) );
    /*    if ( type == IST_PlasticStrainTensor ) {
        answer = status->givePlasDef();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveCumulativePlasticStrain();
        return 1;
    } else if ( ( type == IST_DamageScalar ) || ( type == IST_DamageTensor ) ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
	}*/
}




MisesMatLogStrainStatus :: MisesMatLogStrainStatus(int n, Domain *d, GaussPoint *g) :
  MisesMatStatus(n, d, g)
{
   

    
}

MisesMatLogStrainStatus :: ~MisesMatLogStrainStatus()
{ }

void
MisesMatLogStrainStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
  MisesMatStatus :: printOutputAt(file, tStep);
  fprintf(file, "   Cauchy Stress");
    for ( auto &val : this->CVector ) {
        fprintf(file, "%.4e ", val);
    }
    fprintf(file, "\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void MisesMatLogStrainStatus :: initTempStatus()
{

  MisesMatStatus :: initTempStatus();
  tempCVector      = CVector;
}


// updates internal variables when equilibrium is reached
void
MisesMatLogStrainStatus :: updateYourself(TimeStep *tStep)
{
    MisesMatStatus :: updateYourself(tStep);
    CVector      = tempCVector;
    FVector      = tempFVector;

    
}



} // end namespace oofem
