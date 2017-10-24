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

#ifndef enhancedassumedstrainelementinterface_h
#define enhancedassumedstrainelementinterface_h

#include "interface.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matresponsemode.h"
#include "nlstructuralelement.h"



///@name enhancedassumedstrainelementinterface
//@{
//@}

namespace oofem {

class GaussPoint;
class TimeStep;

/**
 * Enhanced assumed strain element extension interface
 * @author: Martin Horak, nitramkaroh@seznam.cz
 * @reference: Simo, Rifai : A Class of Mixed Assumed Strain Method and The Method of Incompatible Modes
 * @reference: Simo, Armero : Geometrically Non-linear Enahnced Strain Method and The Method of Incompatible Modes
 */

class EnhancedAssumedStrainElementExtensionInterface : public Interface
{



protected:
	FloatArray enahncedDisplacementVector, tempEnahncedDisplacementVector;
	FloatArray internalForcesCorrectionVector, tempInternalForcesCorrectionVector;

	// notation used in reference 1
	FloatMatrix invStiffnessMatrixH, tempInvStiffnessMatrixH;
	FloatMatrix stiffnessMatrixGamma, tempStiffnessMatrixGamma;
	
  
public:
    /**
     * Constructor. Creates element interface belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    EnhancedAssumedStrainElementExtensionInterface(Domain *d);
    /// Destructor.
    virtual ~EnhancedAssumedStrainElementExtensionInterface() { }
	
	void giveEnhancedDisplacementVector(FloatArray &answer){ answer  = this->enahncedDisplacementVector; }
	void giveTempEnhancedDisplacementVector(FloatArray &answer){ answer  = this->tempEnahncedDisplacementVector; }

	void setEnhancedDisplacementVector(FloatArray &eDisplVec){ this->enahncedDisplacementVector = eDisplVec; }
	void setTempEnhancedDisplacementVector(FloatArray &eDisplVec){ this->tempEnahncedDisplacementVector = eDisplVec; }

	void givtInternalForcesCorrectionVector(FloatArray &answer){ answer = internalForcesCorrectionVector; }
	void givTempInternalForcesCorrectionVector(FloatArray &answer){ answer = tempInternalForcesCorrectionVector; }

	void setInternalForcesCorrectionVector(FloatArray &ifcv){ internalForcesCorrectionVector = ifcv; }
	void setTempInternalForcesCorrectionVector(FloatArray &ifcv){ tempInternalForcesCorrectionVector = ifcv; }


	void setInvStiffnessMatrixH(FloatMatrix &smH){ invStiffnessMatrixH = smH; }
	void setTempInvStiffnessMatrixH(FloatMatrix &smH){ tempInvStiffnessMatrixH = smH; }

	void giveInvStiffnessMatrixH(FloatMatrix &answer) { answer = invStiffnessMatrixH; }
	void giveTempInvStiffnessMatrixH(FloatMatrix &answer) { answer = tempInvStiffnessMatrixH; }

	void setStiffnessMatrixGamma(FloatMatrix &smG){ stiffnessMatrixGamma = smG; }
	void setTempStiffnessMatrixGamma(FloatMatrix &smG){ tempStiffnessMatrixGamma = smG; }

	void giveStiffnessMatrixGamma(FloatMatrix &answer){ answer = stiffnessMatrixGamma; }
	void giveTempStiffnessMatrixGamma(FloatMatrix &answer){ answer = tempStiffnessMatrixGamma; }

	void computeEnhancedDisplacementVector(FloatArray &eDispl, StructuralElement *elem, TimeStep* tStep);

	virtual int computeNumberOfEnhancedDofs() = 0;


	void updateYourself();

	/// Enhanced B matrix
	virtual void computeEnhancedBmatrixAt(GaussPoint *gp, FloatMatrix &answer, StructuralElement *elem) = 0;
	/// Enhanced Bh matrix
    virtual void computeEnhancedBHmatrixAt(GaussPoint *gp, FloatMatrix &answer,NLStructuralElement *elem) = 0;
	/// Stiffness correction due to the enhanced strains obtained after static condensation
    void computeStiffnessCorrection(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep, StructuralElement* elem);
	/// Internal forces correction due to the enhanced strains
    void giveInternalForcesCorrectionVector(FloatArray &answer, TimeStep *tStep, StructuralElement* elem);

 protected:   
    /// Right upper block of Stiffness Matrix 
    void computeStiffnessMatrixGamma(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep, StructuralElement* elem);
	/// Right lower block of Stiffness Matrix
    void computeStiffnessMatrixH(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep, StructuralElement* elem);
};

}
#endif
