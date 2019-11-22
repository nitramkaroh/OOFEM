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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "quad1planestrain4eas.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "classfactory.h"
#include <stdio.h>




namespace oofem {

REGISTER_Element( Quad1PlaneStrain4EAS );



Quad1PlaneStrain4EAS :: Quad1PlaneStrain4EAS(int n, Domain *aDomain) :
    Quad1PlaneStrain(n, aDomain),EnhancedAssumedStrainElementExtensionInterface(aDomain)
{
}



void
Quad1PlaneStrain4EAS :: computeEnhancedBmatrixAt(GaussPoint *gp, FloatMatrix &answer,StructuralElement *elem)
//
// Returns the [4x4] enhanced strain - enhanceddisplacement matrix of the receiver,
// evaluated at gp.
{
    
    std::vector<FloatMatrix> E;

    answer.resize(4,4);
    answer.zero();

    

    this-> computeEnhancedModesAt(E, gp);
    
    FloatArray parametricCentroid;
    this->giveElementParametricCentroid(parametricCentroid);
    double j0 = this->giveInterpolation()->giveTransformationJacobian(parametricCentroid, FEIElementGeometryWrapper(this));
    double j = this->giveInterpolation()->giveTransformationJacobian(gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));
    FloatMatrix J0, invJ0;
    this->giveInterpolation()->giveJacobianMatrixAt(J0,parametricCentroid,FEIElementGeometryWrapper(this));
    invJ0.beInverseOf(J0);
    
    FloatMatrix H;

	for (int i = 0; i < E.size(); i++)
	{
		H.beTProductOf(invJ0,E[i]);
		answer.at(1,i+1) += H.at(1,1);
		answer.at(2,i+1) += H.at(2,2);
		answer.at(4,i+1) += H.at(2,1) + H.at(1,2);
		
	}

	answer.times(j0/j);



}


void
Quad1PlaneStrain4EAS :: computeEnhancedBHmatrixAt(GaussPoint *gp, FloatMatrix &answer,NLStructuralElement *elem)
// Returns the [5x8] displacement gradient matrix {BH} of the receiver,
// evaluated at aGaussPoint.
// @todo not checked if correct
{
	  
	

	answer.resize(5,4);
	answer.zero();
	
	std::vector<FloatMatrix> E;
	this-> computeEnhancedModesAt(E, gp);
    
	FloatArray parametricCentroid;
	this->giveElementParametricCentroid(parametricCentroid);
	double j0 = this->giveInterpolation()->giveTransformationJacobian(parametricCentroid, FEIElementGeometryWrapper(this));
	double j = this->giveInterpolation()->giveTransformationJacobian(gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));
	FloatMatrix J0, invJ0;
	this->giveInterpolation()->giveJacobianMatrixAt(J0,parametricCentroid,FEIElementGeometryWrapper(this));
	invJ0.beInverseOf(J0);

	

	for (int i = 0; i < E.size(); i++)
	{
		FloatMatrix H;
		H.beTProductOf(invJ0,E[i]);
		answer.at(1,i+1) += H.at(1,1);
		answer.at(2,i+1) += H.at(2,2);
		answer.at(4,i+1) += H.at(2,1);
		answer.at(5,i+1) += H.at(1,2);
	}

	answer.times(j0/j);
}

	
void 
Quad1PlaneStrain4EAS :: computeEnhancedModesAt(std::vector<FloatMatrix> &answer, GaussPoint* gp)
{
	double ksi, eta;    
	FloatMatrix E;
	
	E.resize(2,2);
	E.zero();

    ksi = (gp->giveNaturalCoordinates()).at(1);
    eta = (gp->giveNaturalCoordinates()).at(2);

	E.at(1,1) = ksi;
	answer.push_back(E);
	
	E.zero();
	E.at(2,2) = eta;
	answer.push_back(E);
	
	E.zero();
	E.at(1,2) = eta;
	answer.push_back(E);

	E.zero();
	E.at(2,1) = ksi;
	answer.push_back(E);
	
	


	

}



} // end namespace oofem
