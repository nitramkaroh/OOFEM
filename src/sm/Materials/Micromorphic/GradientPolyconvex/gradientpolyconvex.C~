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

#include "../sm/Materials/Micromorphic/GradientPolyconvex/gradientpolyconvexmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
  REGISTER_Material(GradientPolyconvexMaterial);

GradientPolyconvexMaterial :: GradientPolyconvexMaterial(int n, Domain *d) :IsotropicLinearElasticMaterial(n, d), MicromorphicMaterialExtensionInterface(d)
{
  Hk = Ak = 0.;
}

GradientPolyconvexMaterial :: ~GradientPolyconvexMaterial()
{ }



void
GradientPolyconvexMaterial :: giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called.");
}



void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix De;
    MaterialMode matMode = gp->giveMaterialMode();
    if (matMode == _PlaneStrain) {
      IsotropicLinearElasticMaterial :: givePlaneStrainStiffMtrx(De, mode, gp, tStep); 
    } else { //@todo check that all other modes are threated as 3d modes
      IsotropicLinearElasticMaterial :: give3dMaterialStiffnessMatrix(De, mode, gp, tStep);
    }
    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
    
    //deformation gradient, its inverse, cofactor, and determinant
    FloatArray vCofF, vInvF;
    FloatMatrix F, invF, cofF;
    F.beMatrixForm(status->giveTempFVector());
    invF.beInverseOf(F);
    double J = F.giveDeterminant();
    cofF.beTranspositionOf(invF);
    cofF.times(J);
    vCofF.beVectorForm(cofF);
    vInvF.beVectorForm(invF);
    // relative stress for cofactor
    FloatArray sCof, mvCof, micromorphicVar;
    micromorphicVar = status->giveTempMicromorphicVar();
    mvCof.beSubArrayOf(micromorphicVar, {2,3,4,5,6,7,8,9,10});
    sCof = vCofF;
    sCof.subtract(mvCof);
    sCof.times(-Hk);
    // relative stress for determinatn
    double sDet = -Hk * (J - micromorphicVar.at(1));

    double gamma = 1;
    
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
	for (int o = 1; o <= 3; o++) {
	  for (int p = 1; p <= 3; p++) {
	    for (int m = 1; m <= 3; m++) {
	      for (int n = 1; n <= 3; n++) {
		answer.at(giveVI(i,j),giveVI(o,p)) += gamma/J * (vInvF.at(giveVI(p,o))*vInvF.at(giveVI(j,i))-vInvF.at(giveVI(j,o))*vInvF.at(giveVI(p,o)))- J*sCof.at(giveVI(m,n))*(vInvF.at(giveVI(j,i))*vInvF.at(giveVI(n,m))-vInvF.at(giveVI(n,i))*vInvF.at(giveVI(j,m)))*vInvF.at(giveVI(p,o))+Hk*J*J*(vInvF.at(giveVI(j,i))*vInvF.at(giveVI(n,m))-vInvF.at(giveVI(n,i))*vInvF.at(giveVI(j,m)))*(vInvF.at(giveVI(p,o))*vInvF.at(giveVI(n,m))-vInvF.at(giveVI(n,o))*vInvF.at(giveVI(p,m)))-J*sCof.at(giveVI(m,n))*(-vInvF.at(giveVI(j,o))*vInvF.at(giveVI(p,i))*vInvF.at(giveVI(n,m))-vInvF.at(giveVI(j,i))*vInvF.at(giveVI(n,o))*vInvF.at(giveVI(p,m))+vInvF.at(giveVI(n,o))*vInvF.at(giveVI(p,i))*vInvF.at(giveVI(j,m))+vInvF.at(giveVI(n,i))*vInvF.at(giveVI(j,o))*vInvF.at(giveVI(p,m)))+(-J*sDet*+Hk*J*J)*vInvF.at(giveVI(j,i))*vInvF.at(giveVI(p,o))-J*sDet*+Hk*J*J*vInvF.at(giveVI(j,o))*vInvF.at(giveVI(p,i));
	      }
	    }
	  }
	}
      }
    }


  
}

void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{


   MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );

   FloatArray vInvF;
   FloatMatrix F, invF;
   F.beMatrixForm(status->giveTempFVector());
   invF.beInverseOf(F);
   double J = F.giveDeterminant();
   vInvF.beVectorForm(invF);

   answer.resize(9,10);
   answer.setColumn(vInvF, 1);
   for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
	for (int o = 1; o <= 3; o++) {
	  for (int p = 1; p <= 3; p++) {
	    answer.at(giveVI(i,j),giveVI(o,p)+1) = (vInvF.at(giveVI(j,i))*vInvF.at(giveVI(p,o))-vInvF.at(giveVI(p,i))*vInvF.at(giveVI(j,o)));
	  }
	}
      }
    }
   answer.times(-J*Hk);   

  
}


void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
  
    FloatArray vInvF;
    FloatMatrix F, invF;
    F.beMatrixForm(status->giveTempFVector());
    invF.beInverseOf(F);
    double J = F.giveDeterminant();
    vInvF.beVectorForm(invF);
    
    answer.resize(9,10);
    answer.setColumn(vInvF, 1);
    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
	for (int o = 1; o <= 3; o++) {
	  for (int p = 1; p <= 3; p++) {
	    answer.at(giveVI(i,j),giveVI(o,p)+1) = (vInvF.at(giveVI(j,i))*vInvF.at(giveVI(p,o))-vInvF.at(giveVI(p,i))*vInvF.at(giveVI(j,o)));
	  }
	}
      }
    }
    answer.times(-J*Hk);   
    
    
}



void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{

    answer.resize(3,3);
    answer.beUnitMatrix();
    answer.times(Hk);

}


void
GradientPolyconvexMaterial :: giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
  MaterialMode matMode = gp->giveMaterialMode();
  if (matMode == _PlaneStrain) {
    answer.resize(6,6);
  } else {
    answer.resize(10,10);
  } 
  answer.beUnitMatrix();
  answer.times(Ak);

}




void
GradientPolyconvexMaterial :: giveFiniteStrainGeneralizedStressVectors(FloatArray &sigma, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &vF, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep)
{
  
    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
    //deformation gradient, its inverse, cofactor, and determinant
    FloatArray vCofF, vInvF;
    FloatMatrix F, invF, cofF;
    F.beMatrixForm(vF);
    invF.beInverseOf(F);
    double J = F.giveDeterminant();
    cofF.beTranspositionOf(invF);
    cofF.times(J);
    vCofF.beVectorForm(cofF);
    vInvF.beVectorForm(invF);
    // relative stress for cofactor
    FloatArray sCof, mvCof;
    mvCof.beSubArrayOf(micromorphicVar, {2,3,4,5,6,7,8,9,10});
    sCof = vCofF;
    sCof.subtract(mvCof);
    sCof.times(-Hk);
    // relative stress for determinatn
    double sDet = -Hk * (J - micromorphicVar.at(1));
    // total relative stress
    s.resize(10);
    s.at(1) = sDet;
    s.addSubVector(sCof, 2);
    // higher order stress
    M = micromorphicVarGrad;
    M.times(Ak);
    


    
    // first PK stress
    FloatArray P(9), c1, c2(9), c3;
    // add Saint-Venant Kirchoff material, i.e. Aijkl Fkl part
    c1 = vCofF;
    double gamma = 1;
    c1.times(-gamma / J / J);

    c3 = vCofF;
    c3.times(-sDet);

    for (int i = 1; i <= 3; i++) {
      for (int j = 1; j <= 3; j++) {
	for (int m = 1; m <=3; m++) {
	  for (int n = 1; n<=3; n++) {
	    c2.at(giveVI(i,j)) = (vInvF.at(giveVI(j,i))*vInvF.at(giveVI(n,m))-vInvF.at(giveVI(n,i))*vInvF.at(giveVI(j,m))) * J * sCof.at(giveVI(m,n));
	  }
	}
      }
    }
    
    P.add(c1);
    P.add(c2);
    P.add(c3);
    


    //    status->letTempStrainVectorBe(displacementGradient);
    status->letTempMicromorphicVarBe(micromorphicVar);
    status->letTempMicromorphicVarGradBe(micromorphicVarGrad);

    status->letTempStressVectorBe(sigma);
    status->letTempMicromorphicStressBe(s);
    status->letTempMicromorphicStressGradBe(M); 
      
}



IRResultType
GradientPolyconvexMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro
    IsotropicLinearElasticMaterial :: initializeFrom(ir);
    
    IR_GIVE_FIELD(ir, Hk, _IFT_MicromorphicMaterialExtensionInterface_Hk);
    IR_GIVE_FIELD(ir, Ak, _IFT_MicromorphicMaterialExtensionInterface_Ak);

    return IRRT_OK;
}

int
GradientPolyconvexMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    MicromorphicMaterialStatus *status = static_cast< MicromorphicMaterialStatus * >( this->giveStatus(gp) );
    

    if( type == IST_MicromorphicStress) {
      StructuralMaterial :: giveFullVectorForm( answer, status->giveStressVector(), gp->giveMaterialMode() );

    } else if( type == IST_MicromorphicStrain ) {
      StructuralMaterial :: giveFullVectorForm( answer, status->giveStrainVector(), gp->giveMaterialMode() );
    } else if( type == IST_MicromorphicRelativeStress ) {
      FloatArray s = status->giveMicromorphicStress();
      // @todo this is correct only for 2d problems, something like giveFullVectorForm for MicromorphicMaterial class would be necessary
      answer.resize(9);
      answer.zero();
      answer.at(6) = s.at(1);
      answer.at(9) = -s.at(1);
    } else if ( type == IST_MicromorphicRelativeStrain ) {
      answer.resize(9);
      answer.zero();
    } else if ( type == IST_MicromorphicHigherOrderStress ) {
      FloatArray M = status->giveMicromorphicStressGrad(); 
      answer.resize(9);
      answer.zero();
      answer.at(6) = M.at(1);
      answer.at(9) = M.at(2);
    } else if ( type == IST_MicromorphicHigherOrderStrain ) {
      FloatArray kappa = status->giveMicromorphicVarGrad();
      answer.resize(9);
      answer.zero();
      answer.at(6) = kappa.at(1);
      answer.at(9) = kappa.at(2);
    } else {
      OOFEM_ERROR("Unknown InternalStateType");
    }
    return 1;
}
    



} // end namespace oofem
