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

#include "../sm/Elements/Beams/beam2beam2danalyticalcontactelement.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "crosssection.h"
#include "node.h"
#include "dof.h"

namespace oofem {
REGISTER_Element(Beam2Beam2dAnalyticalContactElement);


Beam2Beam2dAnalyticalContactElement :: Beam2Beam2dAnalyticalContactElement(int n, Domain *aDomain) : Beam2d(n, aDomain)
{
}
 


void
Beam2Beam2dAnalyticalContactElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    double xc;
    ContactType contact_type;
    this->checkContact(lengthBeam1, lengthBeam2, overlap, xc, tStep, contact_type);  
    if(contact_type == CT_None) {
      answer.clear();
    } else if(contact_type == CT_Beam1_Tip_Beam2) {
      this->computeStiffnessMatrix_TipContact(answer, lengthBeam1, lengthBeam2);
    } else if(contact_type == CT_Beam2_Tip_Beam1) {
      this->computeStiffnessMatrix_TipContact(answer, lengthBeam2, lengthBeam1);
    } else if(contact_type == CT_Overlap) {
      answer.clear();
      //      this->computeStiffnessMatrix_OverlapingContact();
    } else {
      OOFEM_ERROR("Unknown contact type");
    }
           
}



void
Beam2Beam2dAnalyticalContactElement :: computeStiffnessMatrix_TipContact(FloatMatrix &answer, double l1, double l2)
{
      answer = {{0, 1, -l1, 0, -1, -l2},{0, -l1, l1 * l1, 0, l1, l1*l2},{0,  -1, l1, 0, 1, l2},	{0, -l2, l1 * l2, 0, l2, l2*l2}};
      answer.times(3*E*Iy/(l1*l1*l1 + l2*l2*l2) );     
}



void
Beam2Beam2dAnalyticalContactElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{

  double xc;
  ContactType contact_type;
  this->checkContact(lengthBeam1, lengthBeam2, overlap, xc, tStep, contact_type);  
  if(contact_type == CT_None) {
    answer.clear();
  } else if(contact_type == CT_Beam1_Tip_Beam2) {
    this->giveInternalForcesVector_TipContact(answer, lengthBeam1, lengthBeam2, tStep);
  } else if(contact_type == CT_Beam2_Tip_Beam1) {
    this->giveInternalForcesVector_TipContact(answer, lengthBeam2, lengthBeam1, tStep);
  } else if(contact_type == CT_Overlap) {
    this->giveInternalForcesVector_OverlapingContact(answer, lengthBeam1, lengthBeam2, totalLength, tStep);
  } else {
    OOFEM_ERROR("Unknown contact type");
  }

}


void
Beam2Beam2dAnalyticalContactElement :: giveInternalForcesVector_TipContact(FloatArray &answer, double l1, double l2, TimeStep *tStep)
{
    double w1 = this->giveNode(1)->giveDofWithID(2)->giveUnknown(VM_Total, tStep);
    double phi1 = this->giveNode(1)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);
    double w2 = this->giveNode(2)->giveDofWithID(2)->giveUnknown(VM_Total, tStep);
    double phi2 = this->giveNode(2)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);
    
    double F = 3 * E * Iy * (w2 * w1 + phi1 * l1 + phi2 * l2) / ( l1 * l1 *l1 + l2 * l2 * l2);
    answer = {0, -F, F* l1, 0, F, F * l2}; 
}


void
Beam2Beam2dAnalyticalContactElement :: giveInternalForcesVector_OverlapingContact(FloatArray &answer, double l1, double l2, double L, TimeStep *tStep)
{
    double w1 = this->giveNode(1)->giveDofWithID(2)->giveUnknown(VM_Total, tStep);
    double phi1 = this->giveNode(1)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);
    double w2 = this->giveNode(2)->giveDofWithID(2)->giveUnknown(VM_Total, tStep);
    double phi2 = this->giveNode(2)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);
    
    double F = 6 * E * Iy * ( L * ( phi1 + phi2 ) + 2. * (w1 - w2)) / ( L * L * L );
    answer = {0, -F, F * l1, 0, F, F * l2}; 
}



void
Beam2Beam2dAnalyticalContactElement :: checkContact(double l1, double l2, double o, double &xc, TimeStep *tStep, ContactType &contact_type)
{
    double rotationAngle;
    xc = 0;
    this->checkTipContact(rotationAngle, l1, l2 - o, tStep);
    if(rotationAngle <= 0)
      contact_type = CT_Beam1_Tip_Beam2;
    else {
      this->checkTipContact(rotationAngle, l1 - o, l2, tStep);
      if(rotationAngle >= 0){
	contact_type = CT_Beam2_Tip_Beam1;
      } else {
	if(this->checkOverlapingContact(l1,l2, o, xc, tStep)) {
	  contact_type = CT_Overlap;
	} else {
	  contact_type = CT_None;
	}
      }
    }

}

  
void
Beam2Beam2dAnalyticalContactElement :: checkTipContact(double &answer, double l1, double l2, TimeStep *tStep)
{

  double w1 = this->giveNode(1)->giveDofWithID(2)->giveUnknown(VM_Total, tStep);
  double phi1 = this->giveNode(1)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);
  double w2 = this->giveNode(2)->giveDofWithID(2)->giveUnknown(VM_Total, tStep);
  double phi2 = this->giveNode(2)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);
  double F = 3 * E * Iy * (w2 * w1 + phi1 * l1 + phi2 * l2) / ( l1 * l1 *l1 + l2 * l2 * l2);
  answer = F * l1 * l1 / 2 / E / Iy - phi1 - F * l2 * l2 / 2 / E / Iy + phi2;

}



bool
Beam2Beam2dAnalyticalContactElement :: checkOverlapingContact(double l1, double l2, double overlap, double &xc, TimeStep *tStep)
{

  double w1 = this->giveNode(1)->giveDofWithID(2)->giveUnknown(VM_Total, tStep);
  double phi1 = this->giveNode(1)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);
  double w2 = this->giveNode(2)->giveDofWithID(2)->giveUnknown(VM_Total, tStep);
  double phi2 = this->giveNode(2)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);

  double n = 2. * overlap * overlap * phi1 - overlap * l1* phi1 - l1 * l1 * phi1 - 4. * overlap * l2 * phi1 + l1 * l2 * phi1 + 2. * l2 * l2 * phi1 + overlap * overlap * phi2 + overlap * l1 * phi2 - 2. * l1 * l1 * phi2 - 2. * overlap * l2 * phi2 - l1 * l2* phi2 + l2 * l2 * phi2  + 3. * overlap * w1 + 3. * l1* w1 - 3. * l2 * w1 - 3. * overlap * w2 - 3. * l1 * w2 + 3 * l2 * w2;
  double d =  3. * ( overlap - l1 - l2) * ( phi1 + phi2)  + 6. * ( w1- w2);
  xc = n/d;
  if( xc > 0. && xc < overlap)
    return true;
  else {
    return false;
  }
}


      



IRResultType
Beam2Beam2dAnalyticalContactElement :: initializeFrom(InputRecord *ir)
{
    IRResultType result = Beam2d :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    IR_GIVE_FIELD(ir, b1, _IFT_Beam2Beam2dAnalyticalContactElement_b1);
    IR_GIVE_FIELD(ir, b2, _IFT_Beam2Beam2dAnalyticalContactElement_b2);

   
    
    
    return IRRT_OK;
}


void
Beam2Beam2dAnalyticalContactElement ::  postInitialize() 
{
    lengthBeam1  = sqrt(b1.at(1) * b1.at(1) + b1.at(2) * b1.at(2));
    lengthBeam2  = sqrt(b1.at(1) * b1.at(1) + b1.at(2) * b1.at(2));
    
    double dx      = this->giveNode(2)->giveCoordinate(1) - this->giveNode(1)->giveCoordinate(1);
    double dy      = this->giveNode(2)->giveCoordinate(3) - this->giveNode(1)->giveCoordinate(3);
    totalLength =  sqrt(dx * dx + dy * dy);

    //@todo: modify the overlap calculation
    overlap = lengthBeam1 + lengthBeam2 - totalLength;

    FloatArray lc(1);
    Iy   = this->giveCrossSection()->give(CS_InertiaMomentY,  lc, this);
    /*FloatMatrix d;
    this->computeConstitutiveMatrixAt(d, ElasticStiffness, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    E = d.at(1,1);*/
    E = 1;
}


  

} // end namespace oofem
