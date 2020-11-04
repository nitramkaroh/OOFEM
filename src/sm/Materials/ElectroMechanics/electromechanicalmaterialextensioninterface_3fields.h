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

#ifndef electromechanicalmaterialextensioninterface_3field_h
#define electromechanicalmaterialextensioninterface_3field_h

#include "interface.h"
#include "matresponsemode.h"
#include "domain.h"

///@name electromechanicalmaterialextensioninterface
//@{

//@}

namespace oofem {
class FloatMatrix;
class FloatArray;
class GaussPoint;
class TimeStep;



/**
 * Material interface for gradient material models.
 */
class ElectroMechanicalMaterialExtensionInterface : public Interface
{
protected:
    Domain *dom;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    ElectroMechanicalMaterialExtensionInterface(Domain *d){    dom = d;}
    /// Destructor.
    virtual ~ElectroMechanicalMaterialExtensionInterface() { }

  
    virtual void give_FirstPKStressVector_ElectricalDisplacementVector_3d(FloatArray &P, FloatArray &D, GaussPoint *gp, const FloatArray &F, const FloatArray &E, TimeStep *tStep) = 0;

    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) = 0;

    virtual void give3dMaterialStiffnessMatrix_dPdE(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) = 0;
    virtual void give3dMaterialStiffnessMatrix_dDdE(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) = 0;

    
      
};



/**
 * Material interface for gradient material models.
 */
class ElectroMechanicalMaterialExtensionInterface_3Fields : public Interface
{
protected:
    Domain *dom;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    ElectroMechanicalMaterialExtensionInterface_3Fields(Domain *d){    dom = d;}
    /// Destructor.
    virtual ~ElectroMechanicalMaterialExtensionInterface_3Fields() { }

  
    virtual void give_FirstPKStressVector_ElectricFieldVector_3d(FloatArray &P, FloatArray &E, GaussPoint *gp, const FloatArray &F, const FloatArray &D, TimeStep *tStep) = 0;

    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) = 0;

    virtual void give3dMaterialStiffnessMatrix_dPdD(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) = 0;
    virtual void give3dMaterialStiffnessMatrix_dEdD(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep) = 0;

    
      
};
 

}
#endif
