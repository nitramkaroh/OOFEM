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
#include "isohardeningmaterialextensioninterface.h"
#include "inputrecord.h"

namespace oofem {

IsotropicHardeningMaterialExtensionInterface :: IsotropicHardeningMaterialExtensionInterface(Domain *d)  : Interface()
   //
    // constructor
    //
{
  this->domain = d;

}

IRResultType
IsotropicHardeningMaterialExtensionInterface :: initializeFrom(InputRecord *ir)
{

    IRResultType result;              // Required by IR_GIVE_FIELD macro
    // specify the isotropic hardening type
    int isoHardeningTypeRecord = 0; // default
    IR_GIVE_OPTIONAL_FIELD(ir, isoHardeningTypeRecord, _IFT_IsotropicHardeningMaterialExtensionInterface_isohardeningtype);
    if ( isoHardeningTypeRecord == 0 ) {
      this->isoHardeningType = IHT_None;
    } else if ( isoHardeningTypeRecord == 1 ) {
        this->isoHardeningType = IHT_Linear;
    } else if ( isoHardeningTypeRecord == 2 ) {
        this->isoHardeningType = IHT_Saturation;
    } else if ( isoHardeningTypeRecord == 3 ) {
        this->isoHardeningType = IHT_Table;
    } else {
        OOFEM_WARNING("Unknown isohardeningtpye %d", isoHardeningType);
	return IRRT_BAD_FORMAT;
    }
    switch ( isoHardeningType ) {
    case 0:
      break;
    case 1:     
      IR_GIVE_FIELD(ir, H, _IFT_IsotropicHardeningMaterialExtensionInterface_h);
      break;
    case 2:
      IR_GIVE_FIELD(ir, H, _IFT_IsotropicHardeningMaterialExtensionInterface_h);
      IR_GIVE_FIELD(ir, sigInf, _IFT_IsotropicHardeningMaterialExtensionInterface_siginf);
      IR_GIVE_FIELD(ir, delta, _IFT_IsotropicHardeningMaterialExtensionInterface_delta);
      break;
    case 3:
      IR_GIVE_FIELD(ir, kappa,  _IFT_IsotropicHardeningMaterialExtensionInterface_sigh);
      IR_GIVE_FIELD(ir, sigH, _IFT_IsotropicHardeningMaterialExtensionInterface_kappa);
      break;
    default:
      OOFEM_WARNING("Unknown isohardeningtpye %d", isoHardeningType);
      return IRRT_BAD_FORMAT;
    }

    return IRRT_OK;
}
  

double
IsotropicHardeningMaterialExtensionInterface :: giveIsotropicHardeningStress(double kappa)
{

  switch ( isoHardeningType ) {
  case 0:
    return 0;
  case 1:     
    return H*kappa;
    break;
  case 2:
    break;
  case 3:
    break;
  default:
    OOFEM_WARNING("Unknown isohardeningtpye %d", isoHardeningType);
    return 0;
  }
  return 0;
}


double
IsotropicHardeningMaterialExtensionInterface :: giveIsotropicHardeningModulus(double kappa)
{

  switch ( isoHardeningType ) {
  case 0:     
    return 0;
    break;      
  case 1:     
    return H;
    break;
  case 2:
    break;
  case 3:
    break;
  default:
    OOFEM_WARNING("Unknown isohardeningtpye %d", isoHardeningType);
    return IRRT_BAD_FORMAT;
  }
  
  return 0;
}

  

} // end namespace oofem
