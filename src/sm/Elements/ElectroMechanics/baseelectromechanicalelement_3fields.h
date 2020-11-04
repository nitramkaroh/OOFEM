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
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser Base Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser Base Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser Base Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef baseelectromechanicalelement_3fields_h
#define baseelectromechanicalelement_3fields_h

#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/nlstructuralelement.h"

#include "../sm/CrossSections/ElectroMechanics/simpleelectromechanicalcrosssection_3fields.h"

namespace oofem {
/**
 * Base class for electromechanical formulation.
 * @author Martin Horak
 */
  class BaseElectroMechanicalElement_3Fields
{
protected:
    IntArray displacementDofsOrdering, electricDisplacementDofsOrdering, electricPotentialDofsOrdering;
    IntArray locationArray_u, locationArray_phi, locationArray_d;
    

public:
    BaseElectroMechanicalElement_3Fields(int n, Domain *domain);
    virtual ~BaseElectroMechanicalElement_3Fields() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

protected:
  
    /// Pure virtual functions
    virtual NLStructuralElement *giveStructuralElement() = 0;
    virtual void computeElectricPotentialBmatrixAt(GaussPoint *gp, FloatMatrix &Be) = 0;
    virtual void computeElectricDisplacementNmatrixAt(GaussPoint *gp, FloatMatrix &Be) = 0;
    virtual void computeDisplacementFieldBmatrixAt(GaussPoint *gp, FloatMatrix &Bd){this->giveStructuralElement()->computeBmatrixAt(gp, Bd);}

    virtual int giveNumberOfElectricPotentialDofs() = 0;
    virtual int giveNumberOfDisplacementDofs() = 0;
    virtual int giveNumberOfElectricDisplacementDofs() = 0;
    virtual int giveNumberOfDofs() = 0;

    virtual void giveDofManDofIDMask_u(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_phi(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_d(IntArray &answer) = 0;
    /// End of pure virtual functions

    /// @return Reference to the associated crossSection of element.
    SimpleElectroMechanicalCrossSection_3Fields *giveCrossSection(); 

    virtual void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *);
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void compute_FirstPKStressVector_ElectricFieldVector(FloatArray &stress, FloatArray &electricDisplacement, GaussPoint *gp, TimeStep *tStep);
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
    void computeElectricPotentialGradientVector(FloatArray &answer,  GaussPoint *gp, TimeStep *tStep);
    void computeElectricDisplacementVector(FloatArray &answer,  GaussPoint *gp, TimeStep *tStep);   


    void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    void computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);

    virtual IntArray &giveDisplacementDofsOrdering() {return displacementDofsOrdering;}
    virtual IntArray &giveElectricPotentialDofsOrdering() {return electricPotentialDofsOrdering;}
    virtual IntArray &giveElectricDisplacementDofsOrdering() {return electricDisplacementDofsOrdering;}
    void giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_phi, IntArray &locationArray_d, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u,const IntArray &dofIdArray_phi,const IntArray &dofIdArray_d  );
    virtual void postInitialize();
    virtual void updateInternalState(TimeStep *tStep);


};
} // end namespace oofem

#endif
