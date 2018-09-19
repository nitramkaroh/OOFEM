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

#include "ellipticplasticmaterial.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "internalstatetype.h"
#include "contextioerr.h"
#include "intarray.h"
#include "datastream.h"
#include "classfactory.h"
#include "Materials/isolinearelasticmaterial.h"



namespace oofem {
REGISTER_Material(EllipticPlasticMaterial);

EllipticPlasticMaterial :: EllipticPlasticMaterial(int n, Domain *d) : StructuralMaterial(n, d)
{
  linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);

}

EllipticPlasticMaterial :: ~EllipticPlasticMaterial()
{
  delete linearElasticMaterial;
    
}

void
EllipticPlasticMaterial :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseMode mode, GaussPoint *gp,
                                            TimeStep *tStep)
{
    double tempKappa, kappa;
    FloatArray tempTensor2, N;
    EllipticPlasticMaterialStatus *status = static_cast< EllipticPlasticMaterialStatus * >( this->giveStatus(gp) );
    
    FloatArray PS, tempStress, DaN;
    FloatMatrix dN, PP, tempTensor4, Dalgo, cor, De, Ce;
    this->linearElasticMaterial->give3dMaterialStiffnessMatrix(De, TangentStiffness, gp, tStep);
    if ( mode == ElasticStiffness ) {
      answer = De;
    } else if ( mode == TangentStiffness ) {
        kappa = status->giveKappa();
        tempKappa = status->giveTempKappa();
	double dKappa = tempKappa - kappa;
        if ( dKappa > 0. ) {
	  //dN_dSig
	  tempStress = status->giveTempStressVector();
	  this->constructYieldFunctionTensor(PP);
	  // Construction of the derivative of the plastic flow
	  PS.beProductOf(PP, tempStress);
	  FloatArray dFdS;
	  dFdS = PS;
	  dFdS.times(1./sqrt(tempStress.dotProduct(PS)));
	  //////////////////////////////
	  FloatMatrix d2FdS2, temp, tempTensor4;
	  d2FdS2 = PP;
	  temp.beDyadicProductOf(PS, PS);
	  temp.times(1./(tempStress.dotProduct(PS)));
	  d2FdS2.subtract(temp);
	  d2FdS2.times(1./sqrt(tempStress.dotProduct(PS)));
	  // Construction of the gradient Nabla_S of R and SSa tensor
	  Ce.beInverseOf(De);
	  tempTensor4 = d2FdS2;
	  tempTensor4.times(dKappa);
	  tempTensor4.add(Ce);
	  Dalgo.beInverseOf(tempTensor4);
	  double norm = sqrt(tempStress.dotProduct(PS));
	  /// N
	  N = PS;
	  N.times(1./norm);
	  DaN.beProductOf(Dalgo,N);
	  cor.beDyadicProductOf(DaN, DaN);
	  double den = N.dotProduct(DaN);
	  double Hmod = evaluateCurrentPlasticModulus(tempKappa);
	  den += Hmod;
	  cor.times( 1./ den);
	  answer = Dalgo;
	  answer.subtract(cor);	  
        } else {
	  answer = De;
	}
    }

}


double
EllipticPlasticMaterial :: evaluateCurrentYieldStress(GaussPoint *gp)
{
  EllipticPlasticMaterialStatus *status = static_cast< EllipticPlasticMaterialStatus * >( this->giveStatus(gp) );
  return sig0 + H * status->giveTempKappa();
}

double
EllipticPlasticMaterial :: evaluateCurrentPlasticModulus(const double kappa)
{
  return H;
}




void
EllipticPlasticMaterial :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep)
{
    double tempKappa, toSolveScalar;
    FloatArray tempPlasticStrain, tempStress, trialStress, tempTensor2, toSolveTensor, PS, N, incStress;
    FloatMatrix De, Ce, PP, DalgoTensor;

    EllipticPlasticMaterialStatus *status = static_cast< EllipticPlasticMaterialStatus * >( this->giveStatus(gp) );
    FloatMatrix invP(6,6);
    invP.at(1,1) = invP.at(2,2) = invP.at(3,3) = 1;
    invP.at(4,4) = invP.at(5,5) = invP.at(6,6) = 0.5;
    // initialize the plastic strain and cumulative plastic strain
    // by values after the previous step
    tempPlasticStrain = status->givePlasticStrain();
    tempKappa = status->giveKappa();
    // evaluate the trial stress
    this->linearElasticMaterial->give3dMaterialStiffnessMatrix(De, TangentStiffness, gp, tStep);
    Ce.beInverseOf(De);
    trialStress.beProductOf(De, totalStrain - tempPlasticStrain);
    this->constructYieldFunctionTensor(PP);
    // apply the iterative procedure that solves the system of nonlinear equations
    // consisting of the yield condition and discretized flow rule
    // and evaluates:
    // tempKappa   ... cumulative plastic strain at the end of the substep
    // tempStress  ... stress at the end of the substep
    // tempPlasticStrain ... plastic strain at the end of the substep
    tempTensor2.beProductOf(PP, trialStress);
    double SPS = sqrt( trialStress.dotProduct(tempTensor2) );
    double yieldValue = SPS - this->evaluateCurrentYieldStress(gp);
    tempStress = trialStress;
    double rel_yield_tol = 1.e-6;
    if ( yieldValue < rel_yield_tol ) {
      // trial stress in elastic domain
    } else {
        // return to the yield surface needed
        // Initial valuesr
        toSolveTensor.resize(6);
        toSolveTensor.zero();
        toSolveScalar = yieldValue;
        double errorF = yieldValue;
        double errorR = 0;
        DalgoTensor = De;
	PS.beProductOf(PP, tempStress);
	SPS = sqrt(PS.dotProduct(tempStress));
	FloatArray dFdS;
	dFdS = PS;
	dFdS.times(1./sqrt(tempStress.dotProduct(PS)));
        double deltaKappa = 0.;
	bool convergence = false;
	double strain_tol = 1.e-6;
	int flagLoop = 0;
	int max_num_iter = 100;
        do {
  	    double plasModulus = evaluateCurrentPlasticModulus(tempKappa + deltaKappa);
            //*************************************
            //Evaluation of the Recursive Equations
            //*************************************
            tempTensor2.beProductOf(DalgoTensor, dFdS);
            double beta = dFdS.dotProduct(tempTensor2);
            beta +=  plasModulus;
            // Construction of the equation of Delta Kappa
            tempTensor2.beProductOf(DalgoTensor, toSolveTensor);
            double tempScalar = dFdS.dotProduct(tempTensor2);
            double incKappa =  ( toSolveScalar - tempScalar ) / beta;
            tempTensor2 = dFdS;
            tempTensor2.times(incKappa);
            tempTensor2 += toSolveTensor;
            incStress.beProductOf(DalgoTensor, tempTensor2);
	    ////////////////////////////////////////////////////////
	    deltaKappa += incKappa;
	    tempKappa += incKappa;
	    tempStress -= incStress;
	    //*************************************
	    // Evaluation of the f and R
	    //*************************************
	    // Construction of the derivative of the plastic flow
	    PS.beProductOf(PP, tempStress);
	    dFdS = PS;
	    dFdS.times(1./sqrt(tempStress.dotProduct(PS)));
	    //////////////////////////////
	    FloatMatrix d2FdS2, temp, tempTensor4;
	     d2FdS2 = PP;
	     temp.beDyadicProductOf(PS, PS);
	     temp.times(1./(tempStress.dotProduct(PS)));
	     d2FdS2.subtract(temp);
	     d2FdS2.times(1./sqrt(tempStress.dotProduct(PS)));
	    // Construction of the gradient Nabla_S of R and SSa tensor
	    tempTensor4 = d2FdS2;
	    tempTensor4.times(deltaKappa);
	    tempTensor4.add(Ce);
	    DalgoTensor.beInverseOf(tempTensor4);
	    // Evaluation of R
	    tempTensor2.beProductOf( Ce, ( tempStress - trialStress ) );
	    toSolveTensor = tempTensor2 + deltaKappa * dFdS;
	    // Evaluation of f
	    tempTensor2.beProductOf(PP, tempStress);
	    SPS = sqrt( tempStress.dotProduct(tempTensor2) );
	    toSolveScalar = SPS - this->evaluateCurrentYieldStress(gp);
	    //*************************************
	    // Evaluation of the error
	    //*************************************
	    FloatArray errLoop;
	    errLoop = toSolveTensor;
	    errorR = sqrt( errLoop.dotProduct(errLoop) );
	    errorF = fabs(toSolveScalar);
	    
	    flagLoop++;
            convergence = ( fabs(errorF) < rel_yield_tol && errorR < strain_tol );
        } while ( flagLoop <= max_num_iter && !convergence );
    
	status->setTempPlasticStrain(tempPlasticStrain);
	status->setTempKappa(tempKappa);
	status->letTempStressVectorBe(tempStress);
    }
}




void
EllipticPlasticMaterial :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &totalStrain, TimeStep *tStep)
{
    FloatArray effStress, densStress;
    EllipticPlasticMaterialStatus *status = static_cast< EllipticPlasticMaterialStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);
    // compute effective stress using the plasticity model
    performPlasticityReturn(gp, totalStrain, tStep);
    answer =  status->giveTempStressVector();
    // store final damage, strain and stress in status
    status->letTempStrainVectorBe(totalStrain);
}


void
EllipticPlasticMaterial :: constructYieldFunctionTensor(FloatMatrix &answer)
{

  FloatMatrix Id, Iv;
  Id.bePID();
  Id.times( 1.5 * C );
  Iv.beIvol();
  Iv.times( F );
  answer = Id;
  answer.add(Iv);

}





IRResultType
EllipticPlasticMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    result = linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    if ( result != IRRT_OK ) return result;

    sig0 = 0.;
    IR_GIVE_FIELD(ir, sig0, _IFT_EllipticPlasticMaterial_sig0); // uniaxial yield stress

    
    H = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, H, _IFT_EllipticPlasticMaterial_h); // hardening modulus
    C = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, C, _IFT_EllipticPlasticMaterial_c); // hardening modulus
    F = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, F, _IFT_EllipticPlasticMaterial_f); // hardening modulus
    
}



int
EllipticPlasticMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    EllipticPlasticMaterialStatus *status = static_cast< EllipticPlasticMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_PlasticStrainTensor ) {
        answer = status->giveTempPlasDef();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveTempKappa();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}







EllipticPlasticMaterialStatus :: EllipticPlasticMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g)
{
    kappa = tempKappa = 0.0;
    plasticStrain.resize(6);
}



EllipticPlasticMaterialStatus :: ~EllipticPlasticMaterialStatus()
{ }





void
EllipticPlasticMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    /*    fprintf(file, "status { ");
    fprintf( file, "plastrains: %f  %f  %f  %f  %f  %f", this->plasDef.at(1), this->plasDef.at(2), this->plasDef.at(3), this->plasDef.at(4), this->plasDef.at(5), this->plasDef.at(6) );
    fprintf(file, " , kappa %f , dam %f , esed %f , psed %f , tsed %f ", this->tempKappa, this->tempDam, this->tempTSED - this->tempPSED, this->tempPSED, this->tempTSED);
    fprintf(file, "}\n");*/
}


void
EllipticPlasticMaterialStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
    this->tempPlasticStrain = this->plasticStrain;
}


void
EllipticPlasticMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
    this->plasticStrain = this->tempPlasticStrain;
}


contextIOResultType
EllipticPlasticMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    /*    if ( ( iores = plasDef.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream.write(dam) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = effectiveStress.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasFlowDirec.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    */
    return CIO_OK;
}


contextIOResultType
EllipticPlasticMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    /*    // read raw data
    if ( ( iores = plasDef.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream.read(dam) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(beta) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( ( iores = effectiveStress.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = plasFlowDirec.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
    */
    return CIO_OK;
}


MaterialStatus *EllipticPlasticMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new EllipticPlasticMaterialStatus(1, StructuralMaterial :: giveDomain(), gp);
}
} //end namespace oofem
