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

#ifndef incompressiblematerialextensioninterface_h
#define incompressiblematerialextensioninterface_h

#include "interface.h"
#include "matresponsemode.h"



namespace oofem {
class FloatMatrix;
class FloatArray;
class GaussPoint;
class TimeStep;



/**
 * Material interface for treating material incompressibility
 * It is used by mean dilatational method
 * In future, it can be used by mixed formulation using pressure as an independent dof
 */

class IncompressibleMaterialExtensionInterface : public Interface
{
protected:

   
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    IncompressibleMaterialExtensionInterface();
    /// Destructor.
    virtual ~IncompressibleMaterialExtensionInterface() { }
   

    virtual void giveDeviatoricCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep) = 0;
    virtual void giveVolumetricCauchyStressVector_3d(FloatArray &answer, GaussPoint *gp, const double J) = 0;

    virtual void giveVolumetric3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, double J) = 0;    
    virtual void givePressure3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, double J) = 0;
    virtual void giveDeviatoric3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    virtual void giveInitialStiffnessMatrix_Cauchy(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
   

};

}
#endif
