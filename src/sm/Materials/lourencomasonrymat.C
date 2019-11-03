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

#include "lourencomasonrymat.h"
#include "Materials/ortholinearelasticmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include <limits>

namespace oofem {
REGISTER_Material(LourencoMasonryMat);

// constructor
LourencoMasonryMat :: LourencoMasonryMat(int n, Domain *d) : MPlasticMaterial2(n, d)
{
    linearElasticMaterial = new OrthotropicLinearElasticMaterial(n, d);
    this->nsurf = 2;
    // this->rmType = mpm_CuttingPlane;
    this->rmType = mpm_ClosestPoint;

    this->plType = nonassociatedPT; // Rankine-type nonassociated, Hill-type associated
}

// destructor
LourencoMasonryMat :: ~LourencoMasonryMat()
{ }

// specifies whether a given material mode is supported by this model
int
LourencoMasonryMat :: hasMaterialModeCapability(MaterialMode mode)
{
    return mode == _PlaneStress;
}

// reads the model parameters from the input file
IRResultType
LourencoMasonryMat :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    result = StructuralMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    /*
      double value;
    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_ex);
    propertyDictionary.add(Ex, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_ey);
    propertyDictionary.add(Ey, value);

    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_nyxy);
    propertyDictionary.add(NYxy, value);


    IR_GIVE_FIELD(ir, value, _IFT_OrthotropicLinearElasticMaterial_gxy);
    propertyDictionary.add(Gxy, value);
    */


    
    result = linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    if ( result != IRRT_OK ) return result;
    
    // names of orthotropic elastic constants???

    IR_GIVE_FIELD(ir, fcx,_IFT_LourencoMasonryMat_fcx); // Compressive strength in direction x
    IR_GIVE_FIELD(ir, fcy,_IFT_LourencoMasonryMat_fcy); // Compressive strength in direction y
    IR_GIVE_FIELD(ir, ftx,_IFT_LourencoMasonryMat_ftx); // Tensile strength in direction x
    IR_GIVE_FIELD(ir, fty,_IFT_LourencoMasonryMat_fty); // Tensile strength in direction y
    IR_GIVE_FIELD(ir, tauU,_IFT_LourencoMasonryMat_tauU); // Shear strength
    IR_GIVE_FIELD(ir, Gcx,_IFT_LourencoMasonryMat_Gcx); // Fracture energy in compression in direction x
    IR_GIVE_FIELD(ir, Gcy,_IFT_LourencoMasonryMat_Gcy); // Fracture energy in compression in direction y
    IR_GIVE_FIELD(ir, Gtx,_IFT_LourencoMasonryMat_Gtx); // Fracture energy in tension in direction x
    IR_GIVE_FIELD(ir, Gty,_IFT_LourencoMasonryMat_Gty); // Fracture energy in tension in direction y

    beta = -1.00;
    IR_GIVE_OPTIONAL_FIELD(ir, beta,_IFT_LourencoMasonryMat_beta); // Coupling parameter between compressive stresses in two directions
    kappap = 1.2e-3;
    IR_GIVE_OPTIONAL_FIELD(ir, kappap,_IFT_LourencoMasonryMat_kappap); // Value of equivalent plastic strain corresponding to peak of compressive yield stress
    
    return IRRT_OK;
}


MaterialStatus *
LourencoMasonryMat :: CreateStatus(GaussPoint *gp) const
{
    //return new MPlasticMaterial2Status(gp, this->giveSizeOfReducedHardeningVarsVector(gp)); new version
    return new MPlasticMaterial2Status(1, this->giveDomain(), gp, this->giveSizeOfReducedHardeningVarsVector(gp));
    
}

  

// Compute projection array useful for following computations (should I put it in the constructor?)
void
LourencoMasonryMat :: computePiVector(FloatArray &pi) {

  pi.resize( 3 );

  pi.at(1) = 0.5;
  pi.at(2) = 0.5;

  
}

double
LourencoMasonryMat :: computeTensileYieldStress( double E, double f, double G, double kappa, double le ) {

  double answer;

  if ( le > G * E / f / f ) { // edit ft to prevent snap-back at material point level
      f = sqrt( G * E / le );
    }

  answer = f * exp( -1.0 * le * f / G * kappa  );

  return answer;
  
}

double
LourencoMasonryMat :: computeTensileYieldStressKGradient( double E, double f, double G, double kappa, double le ) {

  double answer;

  if ( le > G * E / f / f ) { // edit ft to prevent snap-back at material point level
      f = sqrt( G * E / le );
    }

  answer = -1.0 * le * f * f / G * exp( -1.0 * le * f / G * kappa );

  return answer;
  
}

double
LourencoMasonryMat :: computeCompressiveYieldStress( double E, double f, double G, double kappa, double le ) {

  double answer;
  double sigmap = f;
  double sigmai = f / 3;
  double sigmam = f / 2;
  double sigmar = f / 10;
  double kappam = 75 / 67 * G / le / f + this->kappap;

  if ( kappam < f / E + this->kappap ) { // edit fc to prevent snap-back at material point level
    f = sqrt( 75 / 67 * G * E / le );
  }


  if ( kappa <= this->kappap ) {
    answer  = sigmai + ( sigmap - sigmai ) * sqrt( 2 * kappa / this->kappap - pow( kappa / this->kappap, 2) );
    answer = sigmap / 3. * (1. + 4. * kappa / kappap - 2. * kappa * kappa / kappap / kappap );
    //answer  = sigmai + ( sigmap - sigmai ) / this->kappap * kappa;
  } else if ( kappa <= kappam ) {
    answer = sigmap + ( sigmam - sigmap ) * pow( ( kappa - this->kappap ) / ( kappam - this->kappap ), 2);
    answer = sigmap * (1. - ( (kappa - kappap) / ( kappam - kappap ) ) * ( (kappa - kappap) / ( kappam - kappap ) ) );
  } else {
    answer = sigmar + ( sigmam - sigmar ) * exp( 2 * ( kappa - kappam ) / ( kappam - this->kappap ) * ( sigmam - sigmap ) / ( sigmam - sigmar ) );
    answer = 0;
  }

  return answer;
  
}

double
LourencoMasonryMat :: computeCompressiveYieldStressKGradient( double E, double f, double G, double kappa, double le ) {

  double answer;
  double sigmap = f;
  double sigmai = f / 3;
  double sigmam = f / 2;
  double sigmar = f / 10;
  double kappam = 75 / 67 * G / le / f + this->kappap;

  if ( kappam < f / E + this->kappap ) { // edit fc to prevent snap-back at material point level
    f = sqrt( 75 / 67 * G * E / le );
  }

  if ( kappa <= this->kappap ) {
    //answer = 0.0;//std::numeric_limits<double>::infinity();
    answer = 1.e20;//(sigmap-sigmai)/this->kappap;
    answer = 4. * sigmap / 3. * ( kappap - kappa ) / kappap / kappap;
  } else if ( kappa <= this->kappap ) {
    answer  = ( sigmap - sigmai ) * ( this->kappap - kappa ) / this->kappap /  sqrt( kappa * ( 2 * this->kappap - kappa ));
    answer = -2. * sigmap * ( kappa - kappap ) / ( kappam - kappap )/( kappam - kappap ) ;
  } else if ( kappa <= kappam ) {
    answer = 2 * ( sigmam - sigmap ) * ( kappa - kappap ) / pow( kappam - this->kappap, 2 );
    answer = 0;
  } else {
    answer = 2 * ( sigmam - sigmap ) / ( kappam - this->kappap ) * exp( 2 * ( sigmam - sigmap ) / ( sigmam - sigmar ) * ( kappa - kappam ) / ( kappam - this->kappap ));
    answer = 0;
  }

  return answer;
  
}
 
// Computes backstress (qt or qc) for both yield/loading functions
void
LourencoMasonryMat :: computeTensileBackStressVector(FloatArray &answer, double kappa, double le, GaussPoint *gp ) {

  answer.resize( 3 );
  answer.at(1) = computeTensileYieldStress( linearElasticMaterial->give(Ex, gp), this->ftx, this->Gtx, kappa, le );
  answer.at(2) = computeTensileYieldStress( linearElasticMaterial->give(Ey, gp), this->fty, this->Gty, kappa, le );
   
}

void
LourencoMasonryMat :: computeCompressiveBackStressVector(FloatArray& answer, double kappa, double le, GaussPoint *gp ) {

  answer.resize(2);
  answer.at(1) = computeCompressiveYieldStress( linearElasticMaterial->give(Ex, gp), this->fcx, this->Gcx, kappa, le );
  answer.at(2) = computeCompressiveYieldStress( linearElasticMaterial->give(Ey, gp), this->fcy, this->Gcy, kappa, le );
 
}

void
LourencoMasonryMat :: computeTensileBackStressKGradientVector(FloatArray& answer, double kappa, double le, GaussPoint *gp ) {

   answer.resize(3);
   answer.at(1) = computeTensileYieldStressKGradient( linearElasticMaterial->give(Ex, gp), this->fcx, this->Gcx, kappa, le );
   answer.at(2) = computeTensileYieldStressKGradient( linearElasticMaterial->give(Ey, gp), this->fcy, this->Gcy, kappa, le );
 
}

void
LourencoMasonryMat :: computeCompressiveBackStressKGradientVector(FloatArray& answer, double kappa, double le, GaussPoint *gp ) {

  answer.resize(3);
  answer.at(1) = computeCompressiveYieldStressKGradient( linearElasticMaterial->give(Ex, gp), this->fcx, this->Gcx, kappa, le);
  answer.at(2) = computeCompressiveYieldStressKGradient( linearElasticMaterial->give(Ey, gp), this->fcy, this->Gcy, kappa, le );
 
}
 
// Computes Pt matrix needed for Rankine-type yield function
void 
LourencoMasonryMat :: computePtMatrix(FloatMatrix& answer) {

  double alpha;
  answer.resize(3,3);

  alpha = this->ftx * this->fty / pow( this->tauU, 2);

  answer.at(1,1) = .25;
  answer.at(1,2) = -.25;
  answer.at(2,1) = -.25;
  answer.at(2,2) = .25;
  answer.at(3,3) = alpha;



}

// Computes Pg matrix needed for Rankine-type loading function
void
LourencoMasonryMat :: computePgMatrix(FloatMatrix& answer) {

  answer.resize(3,3);

  answer.at(1,1) = .25;
  answer.at(1,2) = -.25;
  answer.at(2,1) = -.25;
  answer.at(2,2) = .25;
  answer.at(3,3) = 1.0;


}

// Computes Pc matrix needed for Hill-type yield function
void
LourencoMasonryMat :: computePcMatrix(FloatMatrix& answer, double kappa, double le, GaussPoint *gp ) {

  double gamma;
  double sigmacx = computeCompressiveYieldStress( linearElasticMaterial->give(Ex, gp), this->fcx, this->Gcx, kappa, le );
  double sigmacy = computeCompressiveYieldStress( linearElasticMaterial->give(Ey, gp), this->fcy, this->Gcy, kappa, le );
  
  answer.resize( 3, 3 );

  gamma = this->fcx * this->fcy / pow( this->tauU, 2);

  answer.at(1,1) = sigmacy / sigmacx;
  answer.at(1,2) = this->beta / 2;
  answer.at(2,1) = this->beta / 2;
  answer.at(2,2) = sigmacx / sigmacy;
  answer.at(3,3) = gamma;
    
}

// Computes yield function
double
LourencoMasonryMat :: computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables)
{

  // @todo: change to   double le = status->giveLe();
  double le = 1;
  FloatArray pi;
  this->computePiVector(pi);
  FloatArray stressVector;
  StructuralMaterial :: giveReducedSymVectorForm(stressVector, fullStressVector, _PlaneStress);

  
  if ( isurf == 1 ) { //Rankine-type, surface 1
    
    double kappa = strainSpaceHardeningVariables.at( 1 );
    FloatArray backStress;
    this->computeTensileBackStressVector(backStress, kappa, le, gp );
    FloatMatrix Pt;
    this->computePtMatrix(Pt);
    FloatArray eta( 3 );
    FloatArray Pt_eta( 3 );
    
    eta = stressVector - backStress;
    Pt_eta.beProductOf(Pt,eta);
    
    return sqrt(dot( eta, Pt_eta )) + dot( pi, eta );
      
  } else if ( isurf == 2 ) { //Hill-type, surface 2

    double kappa = strainSpaceHardeningVariables.at( 2 );
    FloatArray backStress;
    this->computeCompressiveBackStressVector(backStress, kappa, le, gp);
    FloatMatrix Pc;
    this->computePcMatrix(Pc, kappa, le, gp);
    FloatArray Pc_sigma;
    
    Pc_sigma.beProductOf(Pc,stressVector);
    
    return sqrt(dot( stressVector, Pc_sigma )) - sqrt( backStress.at(1) * backStress.at(2) );
    
  } else {
    
    OOFEM_ERROR("Unknown yield surface %d (only 2 surfaces) ", isurf );
  
  }
}  

//associated and nonassociated flow rule
void
LourencoMasonryMat :: computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &stressSpaceHardeningVars)
{

  switch ( gp->giveMaterialMode() ) {
  case _PlaneStress:
    answer.resize(3);
    break;
  default:
    OOFEM_ERROR("Unknown material mode (%s)", __MaterialModeToString( gp->giveMaterialMode() ) );
  }
  
  answer.zero();
  FloatArray stressVector;
  StructuralMaterial :: giveReducedSymVectorForm(stressVector, fullStressVector, _PlaneStress);
  //@todo: change to double le = status->giveLe(); /// where is giveLe() ?????
  double le = 1;
  if ( isurf == 1 ){ // Rankine-type
    
    FloatMatrix P_t_or_g;

    if ( ftype == 0 ) { // yield function
      this->computePtMatrix(P_t_or_g);
    } else if ( ftype == 1) { // load function
      this->computePgMatrix(P_t_or_g);
    } else {
      OOFEM_ERROR("Unknown function type %d (0 or 1)", ftype );
    }
    // @todo : check this
    //double kappa = strainSpaceHardeningVariables.at( 1 );
    double kappa = stressSpaceHardeningVars.at(1);
    FloatArray backStress;
    this->computeTensileBackStressVector(backStress, kappa, le, gp);
    FloatArray pi;
    this->computePiVector(pi);
    FloatArray eta( 3 );
    
    eta = stressVector - backStress;
    answer.beProductOf( P_t_or_g, eta );
    answer.times( 1. / sqrt(dot( eta, answer ) ) );
    answer.add( pi );
      
  } else if ( isurf == 2 ) {
    //@todo: check this
    //double kappa = strainSpaceHardeningVariables.at( 2 );
    double kappa = stressSpaceHardeningVars.at( 2 );
    FloatArray pi;
    this->computePiVector(pi);
    FloatMatrix Pc;
    this->computePcMatrix(Pc, kappa, le, gp);

    FloatArray Pc_sigma;
    
    Pc_sigma.beProductOf( Pc, stressVector );

    answer = Pc_sigma;
    answer.times( 1. / sqrt( dot( stressVector, Pc_sigma ) ) );
    
  } else {
    
    OOFEM_ERROR("Unknown yield surface %d (only 2 surfaces) ", isurf );
  
  }
}

//necesarry only for mpm_ClosestPoint, see Jirasek: Inelastic analysis of structures, pp. 411.
//Hessian matrix
void
LourencoMasonryMat :: computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables)
{


  FloatArray stressVector;
  StructuralMaterial :: giveReducedSymVectorForm(stressVector, fullStressVector, _PlaneStress);

  switch ( gp->giveMaterialMode() ) {
  case _PlaneStress:
    gradientMatrix.resize(3, 3);
    break;
  default:
    OOFEM_ERROR("Unknown material mode (%s)", __MaterialModeToString( gp->giveMaterialMode() ) );
  }

  gradientMatrix.zero();

  //@todo: change to double le = status->giveLe(); /// where is giveLe() ?????
  double le  = 1;
  if ( isurf == 1) {

    double kappa = strainSpaceHardeningVariables.at( 1 );
    FloatArray backStress;
    computeTensileBackStressVector(backStress, kappa, le, gp);
    FloatMatrix Pg;
    this->computePgMatrix(Pg);
    FloatArray eta( 3 );
    FloatArray Pg_eta( 3 );
    
    eta = stressVector - backStress;
    Pg_eta.beProductOf( Pg, eta );

    gradientMatrix.beDyadicProductOf( Pg_eta, Pg_eta );
    gradientMatrix.times(-1./pow( dot( eta, Pg_eta ), 1.5 ));
    gradientMatrix.add( 1. / sqrt( dot( eta, Pg_eta ) ),  Pg );
      
  } else if ( isurf == 2 ) {

    double kappa = strainSpaceHardeningVariables.at( 2 );
    FloatMatrix Pc;
    this->computePcMatrix(Pc, kappa, le, gp);
    FloatArray Pc_sigma( 3 );

    Pc_sigma.beProductOf( Pc, stressVector );

    gradientMatrix.beDyadicProductOf( Pc_sigma, Pc_sigma );
    gradientMatrix.times( - 1./pow( dot( stressVector, Pc_sigma ), 1.5 ) );
    gradientMatrix.add( 1. / sqrt( dot( stressVector, Pc_sigma ) ), Pc ) ;
    
  } else {
    
    OOFEM_ERROR("Unknown material mode (%s)", __MaterialModeToString( gp->giveMaterialMode() ) );
    
  }
}

void
LourencoMasonryMat :: computeReducedElasticModuli(FloatMatrix &answer,
                                                   GaussPoint *gp,
                                                   TimeStep *tStep)
{  /* Returns elastic moduli in reduced stress-strain space*/
    this->giveLinearElasticMaterial()->giveStiffnessMatrix(answer, ElasticStiffness, gp, tStep);
}

//answer is dkappa (cumulative plastic strain), flow rule
void LourencoMasonryMat :: computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp, const FloatArray &stress, const FloatArray &dlambda, const FloatArray &dplasticStrain, const IntArray &activeConditionMap)
{
    answer.resize(2);
    answer.zero();

    // dkappat = dlambdat
    if ( dlambda.at(1) > 0. ) {
      answer.at(1) = dlambda.at(1);
    }
    
    // dkappac = dlambdac
    if ( dlambda.at(2) > 0. ) {
      answer.at(2) = dlambda.at(2);
    }
}

// Computes the derivative of yield/loading function with respect to kappa_1, kappa_2 etc.
void LourencoMasonryMat :: computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables)
{
    //2 hardening variables for Lourenco masonry model - [kappa_t; kappa_c ]
  
    answer.resize(1);  
    answer.zero();

    FloatArray stressVector;
    StructuralMaterial :: giveReducedSymVectorForm(stressVector, fullStressVector, _PlaneStress);


    //@todo: change to double le = status->giveLe(); /// where is giveLe() ?????
    double le = 1;
    if ( isurf == 1 ) { // Rankine-type surface

      double kappa = strainSpaceHardeningVariables.at( 1 );
      FloatArray backStress;
      this->computeTensileBackStressVector(backStress, kappa, le, gp );
      FloatArray backStressKGradient;
      this->computeTensileBackStressKGradientVector(backStressKGradient, kappa, le, gp );
      FloatMatrix P_t_or_g;
      FloatArray pi;
      this->computePiVector(pi);

      if ( ftype == 0 ) { // yield function
	this->computePtMatrix(P_t_or_g);
      } else if ( ftype == 1) { // load function
	this->computePgMatrix(P_t_or_g);
      } else {
	OOFEM_ERROR("Unknown function type %d (0 or 1)", ftype );
      }

      FloatArray eta = stressVector;
      eta.subtract(backStress);
      
      FloatArray Pt_eta, junk;
      
	
      Pt_eta.beProductOf( P_t_or_g, eta );
      junk = Pt_eta;
      junk.times( 1. / sqrt( dot( eta, Pt_eta ) ) );
      junk.add(pi);

      answer.at(1) = dot( backStressKGradient, junk );
      
    } else if ( isurf == 2 ) {

      double kappa = strainSpaceHardeningVariables.at( 2 );
      FloatArray backStress,backStressKGradient;
      this->computeCompressiveBackStressVector(backStress, kappa, le, gp);
      this->computeCompressiveBackStressKGradientVector(backStressKGradient, kappa, le, gp );
      FloatMatrix Pc;
      this->computePcMatrix(Pc, kappa, le, gp );
      FloatArray Pc_sigma;
      FloatMatrix dPc_dsigmacx(3,3);
      FloatMatrix dPc_dsigmacy(3,3);
      FloatMatrix dPc_dkappa(3,3);
      FloatArray dPc_dkappa_sigma;
      double dfc_dkappac;

      Pc_sigma.beProductOf( Pc, stressVector );

      dPc_dsigmacx.at(1,1) = - backStress.at(2) / pow( backStress.at(1) ,2);
      dPc_dsigmacy.at(1,1) = 1.0 / backStress.at(1);      

      dPc_dsigmacx.at(2,2) = 1.0 / backStress.at(2);
      dPc_dsigmacy.at(2,2) = - backStress.at(1) / pow( backStress.at(2) ,2);  

      dPc_dkappa.add(backStressKGradient.at(1), dPc_dsigmacx);
      dPc_dkappa.add(backStressKGradient.at(2), dPc_dsigmacy);

      dPc_dkappa_sigma.beProductOf( dPc_dkappa, stressVector );

      dfc_dkappac = dot( stressVector, dPc_dkappa_sigma ) / 2 / sqrt(dot( stressVector, Pc_sigma));
      dfc_dkappac -= ( backStress.at(1) * backStressKGradient.at(2) + backStress.at(2) * backStressKGradient.at(1) ) / 2 / sqrt( backStress.at(1) * backStress.at(2) );
            
      answer.at(1) = dfc_dkappac;
      
    } else {
          OOFEM_ERROR("Unknown yield surface %d (only 2 surfaces) ", isurf );
    }
}

//necessary only for mpm_ClosestPoint
//Computes second mixed derivative of loading function with respect to stress and hardening vars.
void LourencoMasonryMat :: computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix, int isurf, GaussPoint *gp, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVariables)
{

  FloatArray stressVector;
  StructuralMaterial :: giveReducedSymVectorForm(stressVector, fullStressVector, _PlaneStress);

  gradientMatrix.resize(3,1); //three stresses (plane stress) and two kappas
  gradientMatrix.zero();

  // @todo: change this to   double le = status->giveLe();
  double le = 1.;

  if ( isurf == 1 ) {

    double kappa = strainSpaceHardeningVariables.at( 1 );
    FloatArray backStress, backStressKGradient;
    this->computeTensileBackStressVector(backStress, kappa, le, gp );
    this->computeTensileBackStressKGradientVector(backStressKGradient, kappa, le, gp );
    FloatArray eta;
    FloatMatrix Pg;
    this->computePgMatrix(Pg);
    FloatMatrix auxMatrix(3,3);
    FloatArray gtSKGradientArray(3);
    FloatArray Pg_eta;
    
    eta = stressVector;
    eta.subtract(backStress);
    Pg_eta.beProductOf( Pg, eta );

    auxMatrix.beDyadicProductOf( Pg_eta, Pg_eta );
    auxMatrix.times(1. / pow( dot( eta, Pg_eta ), 1.5 ) );
    Pg.times( 1. / sqrt( dot( eta, Pg_eta ) ) );
    auxMatrix.subtract( Pg ) ;
    gtSKGradientArray.beProductOf( auxMatrix, backStressKGradient );

    gradientMatrix.at(1,1) = gtSKGradientArray.at(1);
    gradientMatrix.at(2,1) = gtSKGradientArray.at(2);
    gradientMatrix.at(3,1) = gtSKGradientArray.at(3);
    
  } else if ( isurf == 2 ) {

    double kappa = strainSpaceHardeningVariables.at( 2 );
    FloatArray backStress,backStressKGradient;
    this->computeCompressiveBackStressVector(backStress, kappa, le, gp );
    this->computeCompressiveBackStressKGradientVector(backStressKGradient, kappa, le, gp );
    FloatMatrix Pc;
    this->computePcMatrix(Pc, kappa, le, gp );
    FloatArray Pc_sigma;
    FloatMatrix dPc_dsigmacx(3,3);
    FloatMatrix dPc_dsigmacy(3,3);
    FloatMatrix dPc_dkappa(3,3);
    FloatArray dPc_dkappa_sigma;
    FloatArray fcSKGradientArray;

    Pc_sigma.beProductOf( Pc, stressVector );
	    
    dPc_dsigmacx.at(1,1) = - backStress.at(2) / pow( backStress.at(1) ,2);
    dPc_dsigmacx.at(2,2) = 1.0 / backStress.at(2);

    dPc_dsigmacy.at(1,1) = 1.0 / backStress.at(1);
    dPc_dsigmacy.at(2,2) = - backStress.at(1) / pow( backStress.at(2) ,2);

    
    dPc_dkappa.add(backStressKGradient.at(1), dPc_dsigmacx);
    dPc_dkappa.add(backStressKGradient.at(2), dPc_dsigmacy);
    dPc_dkappa_sigma.beProductOf( dPc_dkappa, stressVector );
    
    fcSKGradientArray = dPc_dkappa_sigma;
    fcSKGradientArray.times( 1. / sqrt( dot( stressVector, Pc_sigma ) ) );

    Pc_sigma.times( dot( stressVector, dPc_dkappa_sigma ) / 2 / pow( dot( stressVector, Pc_sigma ), 1.5 ) );
    fcSKGradientArray -= Pc_sigma;


    gradientMatrix.at(1,1) = fcSKGradientArray.at(1);
    gradientMatrix.at(2,1) = fcSKGradientArray.at(2);
    gradientMatrix.at(3,1) = fcSKGradientArray.at(3);
    
  } else {
    OOFEM_ERROR("Unknown yield surface %d (only 2 surfaces) ", isurf );
  }
}

// computes dKappa_i/dsig_j gradient matrix
void LourencoMasonryMat :: computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda)
{
  answer.resize(1,3);
  answer.zero();
}

// computes dKappa_i/dLambda_j for one surface
void LourencoMasonryMat :: computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf, const IntArray &activeConditionMap, const FloatArray &fullStressVector, const FloatArray &strainSpaceHardeningVars, const FloatArray &dlambda)
{

  answer.resize(1,actSurf);
  answer.zero();
  int index = 0;
  for( int i = 1; i <= activeConditionMap.giveSize(); i++) {
    if(activeConditionMap.at(i) == 1) {
      index++;
      answer.at(1, index) = 1;
    }    
  }

  
  /*  for( int i = 1; i <= actSurf; i++) {
    if(activeConditionMap.at(i) == 1) {
      answer.at(1, i) = 1;
    }    
    }*/
  
}

int
LourencoMasonryMat :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    MPlasticMaterial2Status *status = static_cast< MPlasticMaterial2Status * >( this->giveStatus(gp) );
    if ( type == IST_CumPlasticStrain ) {
      FloatArray strainSpaceHardeningVariables;
      strainSpaceHardeningVariables = status->giveStrainSpaceHardeningVarsVector();
      answer.resize(1);
      answer.zero();
      answer.at(1) = strainSpaceHardeningVariables.at(1);
      return 1;
    } else if ( type == IST_CumPlasticStrain_2 ) {
      FloatArray strainSpaceHardeningVariables;
      strainSpaceHardeningVariables = status->giveStrainSpaceHardeningVarsVector();
      answer.resize(1);
      answer.zero();
      answer.at(1) = strainSpaceHardeningVariables.at(2);
      return 1;
    } else {
      return MPlasticMaterial2 :: giveIPValue(answer, gp, type, tStep);
    }
}

  /*
void
LourencoMasonryMat :: giveCharLength(CompoDamageMatStatus *status, GaussPoint *gp, FloatMatrix &elementCs)
{
    FloatArray crackPlaneNormal(2);

    //elementCs.printYourself();

    //normal to x,y,z is the same as in elementCs

    for ( int i = 1; i <= 2; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            crackPlaneNormal.at(j) = elementCs.at(j, i);
        }

        status->elemCharLength.at(i) = gp->giveElement()->giveCharacteristicLength(crackPlaneNormal);
    }
}
  */
 
  
  
} // end namespace oofem
