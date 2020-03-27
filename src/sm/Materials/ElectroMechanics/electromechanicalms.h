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

#ifndef electromechnicalms_h
#define electromechnicalms_h

#include "../sm/Materials/structuralms.h"
#include "floatarray.h"
#include "matstatmapperint.h"

namespace oofem {
class GaussPoint;
class Dictionary;
class Domain;

/**
 * This class implements a structural material status information. It is attribute of
 * gaussPoint. This is only an abstract class, for every instance of material class
 * there should be specialized derived class, which handles are history variables.
 *
 * This is a base class for all material statuses corresponding to materials derived from
 * structural material class.
 * It defines stress and strain vectors and their increments.
 * Functions for accessing these components are defined.
 *
 * Tasks:
 * This is abstract class - only basic functionality is supported like:
 * - maintaining and providing access to stress and strain vectors
 *   (including their increments)
 * - storing and restoring status on tape
 * - printingYourself()
 * - updating Yourself after a new equilibrium state has been reached.
 */
class ElectroMechanicalMaterialStatus : public StructuralMaterialStatus
{
protected:

    /// Equilibrated electrial field vector
    FloatArray EVector;
    /// Temporary electrial field vector(to find balanced state)
    FloatArray tempEVector;
    /// Equilibrated electric displacement vector
    FloatArray DVector;
    /// Temporary electric displacement vector (to find balanced state)
    FloatArray tempDVector;

public:
    /// Constructor. Creates new StructuralMaterialStatus with number n, belonging to domain d and IntegrationPoint g.
    ElectroMechanicalMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~ElectroMechanicalMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    /// Returns the const pointer to receiver's electric field vector
    const FloatArray &giveEVector() const { return EVector; }
    /// Returns the const pointer to receiver's electric displacement
    const FloatArray &giveDVector() const { return DVector; }

    /// Returns the const pointer to receiver's temporary electric field vector
    const FloatArray &giveTempEVector() const { return tempEVector; }
    /// Returns the const pointer to receiver's temporary electric displacement
    const FloatArray &giveTempDVector() const { return tempDVector; }

    
    /// Assigns EVector to given vector v.
    void letEVectorBe(const FloatArray &v) { EVector = v; }
    /// Assigns DVector to given vector v.
    void letDVectorBe(const FloatArray &v) { DVector = v; }
    /// Assigns tempPVector to given vector v
    void letTempEVectorBe(const FloatArray &v) { tempEVector = v; }
    /// Assigns tempFVector to given vector v
    void letTempDVectorBe(const FloatArray &v) { tempDVector = v; }

   
    virtual const char *giveClassName() const { return "ElectroMechanicalMaterialStatus"; }

    /// Functions for MaterialStatusMapperInterface
    virtual void copyStateVariables(const MaterialStatus &iStatus);
    virtual void addStateVariables(const MaterialStatus &iStatus);
};
} // end namespace oofem
#endif // electromechnicalms_h
