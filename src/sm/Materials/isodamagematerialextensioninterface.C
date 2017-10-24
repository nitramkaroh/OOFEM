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


#include "domain.h"
#include "isodamagematerialextensioninterface.h"
#include "inputrecord.h"
#include "mathfem.h"

namespace oofem {

IsotropicDamageMaterialExtensionInterface :: IsotropicDamageMaterialExtensionInterface(Domain *d)  : Interface()
   //
    // constructor
    //
{
  this->domain = d;

}

IRResultType
IsotropicDamageMaterialExtensionInterface :: initializeFrom(InputRecord *ir)
{

    IRResultType result;              // Required by IR_GIVE_FIELD macro
    // specify the isotropic hardening type
    int isoDamageTypeRecord = 0; // default
    IR_GIVE_OPTIONAL_FIELD(ir, isoDamageTypeRecord, _IFT_IsotropicDamageMaterialExtensionInterface_isodamagetype);
    if ( isoDamageTypeRecord == 0 ) {
      this->isoDamageType = IDT_None;
    } else if ( isoDamageTypeRecord == 1 ) {
        this->isoDamageType = IDT_Linear;
    } else if ( isoDamageTypeRecord == 2 ) {
        this->isoDamageType = IDT_Exponential;
    } else if ( isoDamageTypeRecord == 3 ) {
        this->isoDamageType = IDT_Table;
    } else {
        OOFEM_WARNING("Unknown isohardeningtpye %d", isoDamageType);
	return IRRT_BAD_FORMAT;
    }

    switch ( isoDamageType ) {
    case 0:
      break;
    case 1:     
      break;
    case 2:
      IR_GIVE_FIELD(ir, damageCrit, _IFT_IsotropicDamageMaterialExtensionInterface_damageCrit);
      IR_GIVE_FIELD(ir, damageExp,  _IFT_IsotropicDamageMaterialExtensionInterface_damageExp);
      break;
    case 3:
      IR_GIVE_FIELD(ir, kappaTab,  _IFT_IsotropicDamageMaterialExtensionInterface_kappaTab);
      IR_GIVE_FIELD(ir, damageTab, _IFT_IsotropicDamageMaterialExtensionInterface_damageTab);
      break;
    default:
      OOFEM_WARNING("Unknown isohardeningtpye %d", isoDamageType);
      return IRRT_BAD_FORMAT;
    }

    return IRRT_OK;
}



double
IsotropicDamageMaterialExtensionInterface :: giveDamageParam(double tempKappa)
{

  if ( tempKappa > 0. ) {
    switch ( isoDamageType ) {
    case 0:
      return 0;
    case 1:     
      return 0;
      break;
    case 2:
      return damageCrit * ( 1.0 - exp(-damageExp * tempKappa) );
      break;
    case 3:
      return 0;
      break;
    default:
      OOFEM_WARNING("Unknown isohardeningtpye %d", isoDamageType);
      return 0;
    } 
  } else {
    return 0.;
  }
}

double
IsotropicDamageMaterialExtensionInterface :: giveDamageParamPrime(double tempKappa)
{

  if ( tempKappa > 0. ) {
    switch ( isoDamageType ) {
    case 0:
      return 0;
    case 1:     
      return 0;
      break;
    case 2:
      return damageCrit * damageExp * exp(-damageExp * tempKappa);
      break;
    case 3:
      return 0;
      break;
    default:
      OOFEM_WARNING("Unknown isohardeningtpye %d", isoDamageType);
      return 0;
    } 
  } else {
    return 0.;
  }
  
}

  

  

} // end namespace oofem
