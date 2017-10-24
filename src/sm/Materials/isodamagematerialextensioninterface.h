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
 *               Copyright (C) 1993 - 2017   Borek Patzak
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

#ifndef isodamagematerialextensioninterface_h
#define isodamagematerialextensioninterface_h

#include "interface.h"
#include "matresponsemode.h"
#include "floatarray.h"
#include "inputrecord.h"


///@name graddpdmaterialextensioninterface
//@{
#define _IFT_IsotropicDamageMaterialExtensionInterface_isodamagetype "isodamagetype"

#define _IFT_IsotropicDamageMaterialExtensionInterface_damageCrit "damagecrit"
#define _IFT_IsotropicDamageMaterialExtensionInterface_damageExp "damageexp"

#define _IFT_IsotropicDamageMaterialExtensionInterface_kappaTab "kappatab"
#define _IFT_IsotropicDamageMaterialExtensionInterface_damageTab "damagetab"


//@}

namespace oofem {
class FloatMatrix;
class GaussPoint;
class TimeStep;

/**
 * Material interface for gradient material models.
 */
class IsotropicDamageMaterialExtensionInterface : public Interface
{
protected:
    Domain *domain;

    /** Variable characterizing the isotropic hardening type
     */
   enum IsotropicDamageType {
     IDT_None=0,
     IDT_Linear=1,
     IDT_Exponential=2,
     IDT_Table=3,
     IDT_Unknown = 100
   };
    
   IsotropicDamageType isoDamageType;
   /**
     * Maximum damage that can be reached 
     */
    double damageCrit;

    /**
     * Damage exponent
     */
    double damageExp;

    /**
     * Array of values of accumulated plastic strain
     * Used when the damage function is given by table values
     */
    FloatArray kappaTab;

    
    /**
     * Array of values of damage
     * Used when the damage function is given by table values
     */
    FloatArray damageTab;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    IsotropicDamageMaterialExtensionInterface(Domain *d);
    /// Destructor.
    virtual ~IsotropicDamageMaterialExtensionInterface(){}
    virtual double giveDamageParam(double kappa);
    virtual double giveDamageParamPrime(double kappa);
    
    virtual IRResultType initializeFrom(InputRecord *ir);


};


}
#endif
