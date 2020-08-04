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

#include "../sm/CrossSections/ElectroMechanics/simpleelectromechanicalcrosssection_3fields.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Materials/ElectroMechanics/electromechanicalmaterialextensioninterface_3fields.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "engngm.h"

namespace oofem {
REGISTER_CrossSection(SimpleElectroMechanicalCrossSection);


void
SimpleElectroMechanicalCrossSection_3Fields :: give_FirstPKStressVector_ElectricFieldVector(FloatArray &P, FloatArray &E, GaussPoint *gp, const FloatArray &F, const FloatArray &D, TimeStep *tStep)
{
    // This function returns the first Piola-Kirchoff stress in vector format and vector of electrical displacement
    // corresponding to a given deformation gradient according to the stress-deformation
    // mode stored in the each gp.

    MaterialMode mode = gp->giveMaterialMode();
    ElectroMechanicalMaterialExtensionInterface_3Fields *elmechMat = static_cast< ElectroMechanicalMaterialExtensionInterface_3Fields * >(this->giveMaterialInterface(ElectroMechanicalMaterialExtensionInterface_3FieldsType, gp) );
      if ( !elmechMat ) {
        OOFEM_ERROR("Material doesn't implement the required Electro-Mechanical Material interface!");
      }

    
    if ( mode == _3dMat ) {
      elmechMat->give_FirstPKStressVector_ElectricFieldVector_3d(P, E, gp, F, D, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}




void
SimpleElectroMechanicalCrossSection_3Fields :: giveElectroMechanicalConstitutiveMatrix_dPdF(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    ElectroMechanicalMaterialExtensionInterface_3Fields *elmechMat = static_cast< ElectroMechanicalMaterialExtensionInterface_3Fields * >(this->giveMaterialInterface(ElectroMechanicalMaterialExtensionInterface_3FieldsType, gp) );
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        elmechMat->give3dMaterialStiffnessMatrix_dPdF(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
SimpleElectroMechanicalCrossSection_3Fields :: giveElectroMechanicalConstitutiveMatrix_dPdD(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    ElectroMechanicalMaterialExtensionInterface_3Fields *elmechMat = static_cast< ElectroMechanicalMaterialExtensionInterface_3Fields * >(this->giveMaterialInterface(ElectroMechanicalMaterialExtensionInterfaceType, gp) );

    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        elmechMat->give3dMaterialStiffnessMatrix_dPdD(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}

  
void
SimpleElectroMechanicalCrossSection_3Fields :: giveElectroMechanicalConstitutiveMatrix_dEdD(FloatMatrix &answer,
                                               MatResponseMode rMode, GaussPoint *gp,
                                               TimeStep *tStep)
{
    ElectroMechanicalMaterialExtensionInterface_3Fields *elmechMat = static_cast< ElectroMechanicalMaterialExtensionInterface_3Fields * >(this->giveMaterialInterface(ElectroMechanicalMaterialExtensionInterface_3FieldsType, gp) );

    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        elmechMat->give3dMaterialStiffnessMatrix_dEdD(answer, rMode, gp, tStep);
    } else {
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mode) );
    }
}

  

IRResultType
SimpleElectroMechanicalCrossSection_3Fields :: initializeFrom(InputRecord *ir)
//
// instanciates receiver from input record
//
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    double thick = 0.0;
    if ( ir->hasField(_IFT_SimpleElectroMechanicalCrossSection_3Fields_thick) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, thick, _IFT_SimpleElectroMechanicalCrossSection_3Fields_thick);
        propertyDictionary.add(CS_Thickness, thick);
    }

    double width = 0.0;
    if ( ir->hasField(_IFT_SimpleElectroMechanicalCrossSection_3Fields_width) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, width, _IFT_SimpleElectroMechanicalCrossSection_3Fields_width);
        propertyDictionary.add(CS_Width, width);
    }

    double area = 0.0;
    if ( ir->hasField(_IFT_SimpleElectroMechanicalCrossSection_3Fields_area) ) {
        IR_GIVE_FIELD(ir, area, _IFT_SimpleElectroMechanicalCrossSection_3Fields_area);
    } else {
        area = thick * width;
    }
    propertyDictionary.add(CS_Area, area);


    this->materialNumber = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->materialNumber, _IFT_SimpleElectroMechanicalCrossSection_3Fields_MaterialNumber);
  
    return CrossSection :: initializeFrom(ir);
}


void
SimpleElectroMechanicalCrossSection_3Fields :: createMaterialStatus(GaussPoint &iGP)
{
    Material *mat = domain->giveMaterial(materialNumber);
    MaterialStatus *matStat = mat->CreateStatus(& iGP);
    iGP.setMaterialStatus(matStat);
}


bool
SimpleElectroMechanicalCrossSection_3Fields :: isCharacteristicMtrxSymmetric(MatResponseMode rMode)
{
    if ( this->giveMaterialNumber() ) {
        return this->domain->giveMaterial( this->giveMaterialNumber() )->isCharacteristicMtrxSymmetric(rMode);
    } else {
        return false; // Bet false...
    }
}


Material
*SimpleElectroMechanicalCrossSection_3Fields :: giveMaterial(IntegrationPoint *ip)
{
    if ( this->giveMaterialNumber() ) {
        return this->giveDomain()->giveMaterial( this->giveMaterialNumber() );
    } else {
        return ip->giveElement()->giveMaterial();
    }
}


double
SimpleElectroMechanicalCrossSection_3Fields :: give(int aProperty, GaussPoint *gp)
{
    return this->giveMaterial(gp)->give(aProperty, gp);
}


int
SimpleElectroMechanicalCrossSection_3Fields :: giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_CrossSectionNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    }
    return this->giveMaterial(ip)->giveIPValue(answer, ip, type, tStep);
}



int
SimpleElectroMechanicalCrossSection_3Fields :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    Material *mat = this->giveDomain()->giveMaterial(this->materialNumber);
    
    if ( !dynamic_cast< StructuralMaterial * >(mat) ) {
        OOFEM_WARNING( "material %s without structural support", mat->giveClassName() );
        result = 0;
    }

    return result;
}


void
SimpleElectroMechanicalCrossSection_3Fields :: giveInputRecord(DynamicInputRecord &input)
{
    CrossSection :: giveInputRecord(input);

    if ( this->propertyDictionary.includes(CS_Thickness) ) {
        input.setField(this->propertyDictionary.at(CS_Thickness), _IFT_SimpleElectroMechanicalCrossSection_3Fields_thick);
    }

    if ( this->propertyDictionary.includes(CS_Width) ) {
        input.setField(this->propertyDictionary.at(CS_Width), _IFT_SimpleElectroMechanicalCrossSection_3Fields_width);
    }

    if ( this->propertyDictionary.includes(CS_Area) ) {
        input.setField(this->propertyDictionary.at(CS_Area), _IFT_SimpleElectroMechanicalCrossSection_3Fields_area);
    }

    
    input.setField(this->materialNumber, _IFT_SimpleElectroMechanicalCrossSection_3Fields_MaterialNumber);
    
}


  
Interface
*SimpleElectroMechanicalCrossSection_3Fields :: giveMaterialInterface(InterfaceType t, IntegrationPoint *ip)
{
    return this->giveMaterial(ip)->giveInterface(t);
}



int
SimpleElectroMechanicalCrossSection_3Fields :: packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->packUnknowns(buff, tStep, gp);
}

int
SimpleElectroMechanicalCrossSection_3Fields :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)
{
    return this->giveMaterial(gp)->unpackAndUpdateUnknowns(buff, tStep, gp);
}

int
SimpleElectroMechanicalCrossSection_3Fields :: estimatePackSize(DataStream &buff, GaussPoint *gp)
{
    return this->giveMaterial(gp)->estimatePackSize(buff, gp);
}
} // end namespace oofem
