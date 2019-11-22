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
    if(E == -1) {
      E = this->initYoungModulus(tStep);
    }
    
    if(contact_type == CT_None) {
      answer.clear();
    } else if(contact_type == CT_Beam1_Tip_Beam2) {
      this->computeStiffnessMatrix_TipContact(answer, lengthBeam1, lengthBeam2 - overlap);
    } else if(contact_type == CT_Beam2_Tip_Beam1) {
      this->computeStiffnessMatrix_TipContact(answer, lengthBeam1 - overlap, lengthBeam2);
    } else if(contact_type == CT_Overlap) {
      answer.clear();
      this->computeStiffnessMatrix_OverlapingContact(answer, totalLength);
    } else {
      OOFEM_ERROR("Unknown contact type");
    }
           
}



void
Beam2Beam2dAnalyticalContactElement :: computeStiffnessMatrix_TipContact(FloatMatrix &answer, double l1, double l2)
{
  answer = {{0, 0, 0, 0, 0, 0}, {0, 1, -l1, 0, -1, -l2},{0, -l1, l1 * l1, 0, l1, l1*l2},{0, 0, 0, 0, 0, 0}, {0,  -1, l1, 0, 1, l2}, {0, -l2, l1 * l2, 0, l2, l2*l2}};
      answer.times(3*E*Iy/(l1*l1*l1 + l2*l2*l2) );     
}

void
Beam2Beam2dAnalyticalContactElement :: computeStiffnessMatrix_OverlapingContact(FloatMatrix &answer, double L)
{
  
  double EI = E*Iy;
  double L3 = L * L * L;
  double L2 = L * L;
  
  answer = {{0, 0, 0, 0, 0, 0}, {0, 12*EI/L3, -6*EI/L2, 0, -12*EI/L3, -6*EI/L2},{0, -6*EI/L2, 4*EI/L, 0, 6*EI/L2, 2*EI/L},{0, 0, 0, 0, 0, 0}, {0, -12*EI/L3, 6*EI/L2, 0, 12*EI/L3, 6*EI/L2}, {0, -6*EI/L2, 2*EI/L, 0, 6*EI/L2, 4*EI/L}};
}



void
Beam2Beam2dAnalyticalContactElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
  if(E == -1) {
      E = this->initYoungModulus(tStep);
    }
  double xc;
  ContactType contact_type;
  this->checkContact(lengthBeam1, lengthBeam2, overlap, xc, tStep, contact_type);  
  if(contact_type == CT_None) {
    answer.clear();
  } else if(contact_type == CT_Beam1_Tip_Beam2) {
    this->giveInternalForcesVector_Contact(answer, lengthBeam1, lengthBeam2 - overlap, tStep);
  } else if(contact_type == CT_Beam2_Tip_Beam1) {
    this->giveInternalForcesVector_Contact(answer, lengthBeam1 - overlap, lengthBeam2, tStep);
  } else if(contact_type == CT_Overlap) {
    this->giveInternalForcesVector_Contact(answer, lengthBeam1 - xc, lengthBeam2 - overlap + xc , tStep);
  } else {
    OOFEM_ERROR("Unknown contact type");
  }

}


void
Beam2Beam2dAnalyticalContactElement :: giveInternalForcesVector_Contact(FloatArray &answer, double l1, double l2, TimeStep *tStep)
{
    // compute position after deformation
    double w1, phi1, w2, phi2;
    this->giveLocalDofs(w1, phi1, w2, phi2, tStep);
    
    double F = 3 * E * Iy * (w2 - w1 + phi1 * l1 + phi2 * l2) / ( l1 * l1 *l1 + l2 * l2 * l2);
    answer = {0, -F, F* l1, 0, F, F * l2}; 
}



void
Beam2Beam2dAnalyticalContactElement :: checkContact(double l1, double l2, double o, double &xc, TimeStep *tStep, ContactType &contact_type)
{
    double rotationAngle;
    FloatArray x, x1, x2;
    // checking beam positions
    x1 = *this->giveNode(1)->giveCoordinates();
    x2 = *this->giveNode(2)->giveCoordinates();
    x = x2 + b2;
    double sign = -n1.dotProduct(x-x1);
    if(sign == 0) {
      sign  = 1;
    } else {
      sign /= fabs(sign);
    }
    // compute position after deformation
    double w1, phi1, w2, phi2;
    this->giveLocalDofs(w1, phi1, w2, phi2, tStep);
    // check two contact positions
    double cond1 = sign * ( phi1 * lengthBeam1 + phi2 * ( lengthBeam2 - overlap ) + w2 - w1 );
    double cond2 = sign * ( phi2 * lengthBeam2 + phi1 * ( lengthBeam1 - overlap ) + w2 - w1 );
    if (cond1 <= 0 && cond2 <= 0) {
      contact_type = CT_None;
    } else {
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
	  }
	}
      }
    }
}



  
  
void
Beam2Beam2dAnalyticalContactElement :: checkTipContact(double &answer, double l1, double l2, TimeStep *tStep)
{

  double w1, phi1, w2, phi2;
  this->giveLocalDofs(w1, phi1, w2, phi2, tStep);
  double F = 3 * E * Iy * (w2 - w1 + phi1 * l1 + phi2 * l2) / ( l1 * l1 *l1 + l2 * l2 * l2);
  answer = F * l1 * l1 / 2 / E / Iy - phi1 - F * l2 * l2 / 2 / E / Iy + phi2;

}



bool
Beam2Beam2dAnalyticalContactElement :: checkOverlapingContact(double l1, double l2, double overlap, double &xc, TimeStep *tStep)
{

  double w1 = this->giveNode(1)->giveDofWithID(3)->giveUnknown(VM_Total, tStep);
  double phi1 = this->giveNode(1)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);
  double w2 = this->giveNode(2)->giveDofWithID(3)->giveUnknown(VM_Total, tStep);
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


      
void
Beam2Beam2dAnalyticalContactElement :: giveLocalDofs(double &w1, double &phi1, double &w2, double &phi2, TimeStep *tStep)
{

  FloatArray global_Dofs(6), local_Dofs;
  global_Dofs.at(1) = this->giveNode(1)->giveDofWithID(1)->giveUnknown(VM_Total, tStep);
  global_Dofs.at(2) = this->giveNode(1)->giveDofWithID(3)->giveUnknown(VM_Total, tStep);
  global_Dofs.at(3) = this->giveNode(1)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);
  global_Dofs.at(4) = this->giveNode(2)->giveDofWithID(1)->giveUnknown(VM_Total, tStep);
  global_Dofs.at(5) = this->giveNode(2)->giveDofWithID(3)->giveUnknown(VM_Total, tStep);
  global_Dofs.at(6) = this->giveNode(2)->giveDofWithID(5)->giveUnknown(VM_Total, tStep);

  FloatMatrix R;
  this->computeGtoLRotationMatrix(R);
  local_Dofs.beProductOf(R, global_Dofs);
  w1 = local_Dofs.at(2);
  phi1 = local_Dofs.at(3);
  w2 = local_Dofs.at(5);
  phi2 = local_Dofs.at(6);

    
}

bool
Beam2Beam2dAnalyticalContactElement :: computeGtoLRotationMatrix(FloatMatrix &answer)
{

   answer.resize(6, 6);
   answer.zero();
  
   answer.at(1, 1) =  b1.at(1) / lengthBeam1;
   answer.at(1, 2) =  b1.at(3) / lengthBeam1;
   answer.at(2, 1) =  -b1.at(3) / lengthBeam1;
   answer.at(2, 2) =  b1.at(1) / lengthBeam1;
   answer.at(3, 3) =  1.;
   answer.at(4, 4) = -b2.at(1) / lengthBeam2;
   answer.at(4, 5) = -b2.at(3) / lengthBeam2;
   answer.at(5, 4) =  b2.at(3) / lengthBeam2;
   answer.at(5, 5) = -b2.at(1) / lengthBeam2;
   answer.at(6, 6) =  1.;
   return true;
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

double
Beam2Beam2dAnalyticalContactElement :: initYoungModulus( TimeStep *tStep)
{

    FloatMatrix d;
    GaussPoint *gp = NULL;
    this->giveStructuralCrossSection()->giveStiffnessMatrix_1d(d, ElasticStiffness, gp, tStep);
    return d.at(1,1);
}



void
Beam2Beam2dAnalyticalContactElement ::  postInitialize() 
{
    lengthBeam1  = sqrt(b1.at(1) * b1.at(1) + b1.at(3) * b1.at(3));
    lengthBeam2  = sqrt(b2.at(1) * b2.at(1) + b2.at(3) * b2.at(3));

    FloatArray x1, x2, dx;
    x1 = *this->giveNode(1)->giveCoordinates();
    x2 = *this->giveNode(2)->giveCoordinates();
    dx = x2 - x1;
    totalLength =  sqrt(dx.at(1) * dx.at(1) + dx.at(3) * dx.at(3));

    n1 = { -b1.at(3) / lengthBeam1, 0, b1.at(1) / lengthBeam1 };
    n2 = { -b2.at(3) / lengthBeam2, 0, b2.at(1) / lengthBeam2 };

    //@todo: check this
    overlap = b1.dotProduct(x1 + b1 - x2 - b2)/lengthBeam1;
   
    FloatArray lc(1);
    Iy   = this->giveCrossSection()->give(CS_InertiaMomentY,  lc, this);

}



  void
Beam2Beam2dAnalyticalContactElement :: printOutputAt(FILE *File, TimeStep *tStep)
{
  Beam2d :: printOutputAt(File, tStep);
  ContactType contact_type;
  double xc;
  this->checkContact(lengthBeam1, lengthBeam2, overlap, xc, tStep, contact_type);  
  if(contact_type == CT_None) {
    fprintf(File, "no contact\n");
  } else if (contact_type == CT_Beam1_Tip_Beam2) {
    fprintf(File, "tip of beam1 on beam2 contact\n");
  }  else if (contact_type == CT_Beam2_Tip_Beam1) {
    fprintf(File, "tip of beam2 on beam1 contact\n");
  } else if (contact_type == CT_Overlap) {
    fprintf(File, "overlapping contact\n");
  }
    
}
  
  

} // end namespace oofem
