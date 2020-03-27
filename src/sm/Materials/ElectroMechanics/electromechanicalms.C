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
#include "../sm/Materials/ElectroMechanics/electromechanicalms.h"
#include "contextioerr.h"
#include "gausspoint.h"

namespace oofem {
  
  ElectroMechanicalMaterialStatus :: ElectroMechanicalMaterialStatus(int n, Domain *d, GaussPoint *g) :   StructuralMaterialStatus(n, d, g),  EVector(), tempEVector(),  DVector(), tempDVector()
{
    /// Hack to prevent crashing when there are only "loose" gausspoints
    if ( gp->giveIntegrationRule() == NULL ) {
        return;
    }
    EVector.resize(3);
    DVector.resize(3);
    
    tempEVector = EVector;
    tempDVector = DVector;

}


ElectroMechanicalMaterialStatus :: ~ElectroMechanicalMaterialStatus() { }


void ElectroMechanicalMaterialStatus :: printOutputAt(FILE *File, TimeStep *tStep)
// Prints the strains and stresses on the data file.
{
    StructuralMaterialStatus :: printOutputAt(File, tStep);
    fprintf(File, "\n");
}


void ElectroMechanicalMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    StructuralMaterialStatus :: updateYourself(tStep);
    EVector      = tempEVector;
    DVector      = tempDVector;
}


void ElectroMechanicalMaterialStatus :: initTempStatus()
//
// initialize record at the begining of new load step
//
{
    StructuralMaterialStatus :: initTempStatus();

    tempEVector      = EVector;
    tempDVector      = DVector;
    
}


contextIOResultType
ElectroMechanicalMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
{
    contextIOResultType iores;

    return CIO_OK;
}


contextIOResultType
ElectroMechanicalMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    return CIO_OK;
}

void ElectroMechanicalMaterialStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const ElectroMechanicalMaterialStatus &status = dynamic_cast< ElectroMechanicalMaterialStatus & >(tmpStat);

    EVector = status.giveEVector();
    tempEVector = status.giveTempEVector();
    DVector = status.giveDVector();
    tempDVector = status.giveTempDVector();
}

void ElectroMechanicalMaterialStatus :: addStateVariables(const MaterialStatus &iStatus)
{
    printf("Entering ElectroMechanicalMaterialStatus :: addStateVariables().\n");
}

  
} // end namespace oofem

