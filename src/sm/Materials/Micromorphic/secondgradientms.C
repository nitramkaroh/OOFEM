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

#include "../sm/Materials/Micromorphic/secondgradientms.h"
#include "../sm/Materials/structuralmaterial.h"
#include "contextioerr.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "gausspoint.h"

namespace oofem {
  SecondGradientMaterialStatus :: SecondGradientMaterialStatus(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g), micromorphicVar(), micromorphicVarGrad(), micromorphicStressGrad(), tempMicromorphicVar(), tempMicromorphicVarGrad(), tempMicromorphicStressGrad() 
{
    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    strainVector.resize(rsize);
    stressVector.resize(rsize);
    // reset temp vars.
    tempStressVector = stressVector;
    tempStrainVector = strainVector;
    /// initialization of micromorphic vars

    micromorphicVar.resize(1);
    micromorphicVarGrad.resize(2);

    tempMicromorphicVar = micromorphicVar;
    tempMicromorphicVarGrad = micromorphicVarGrad;


    micromorphicStressGrad.resize(2);
    tempMicromorphicStressGrad = micromorphicStressGrad;

     
}


SecondGradientMaterialStatus :: ~SecondGradientMaterialStatus() { }


void SecondGradientMaterialStatus :: printOutputAt(FILE *File, TimeStep *tStep)
// Prints the strains and stresses on the data file.
{
    FloatArray helpVec;
    MaterialStatus :: printOutputAt(File, tStep);

    fprintf(File, "  strains ");
    // here should be sym for microdil and without sym for cosserat!!
    //    StructuralMaterial :: giveFullSymVectorForm( helpVec, strainVector, gp->giveMaterialMode() );
    StructuralMaterial :: giveFullVectorForm( helpVec, strainVector, gp->giveMaterialMode() );
    for ( auto &var : helpVec ) {
      fprintf( File, " %.4e", var );
    }
    
    fprintf(File, "\n              stresses");
    //StructuralMaterial :: giveFullSymVectorForm( helpVec, stressVector, gp->giveMaterialMode() );

    StructuralMaterial :: giveFullVectorForm( helpVec, stressVector, gp->giveMaterialMode() );
    
    for ( auto &var : helpVec ) {
      fprintf( File, " %.4e", var );
    }
    fprintf(File, "\n");

}




void SecondGradientMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    StructuralMaterialStatus :: updateYourself(tStep);

    micromorphicVar = tempMicromorphicVar;
    micromorphicVarGrad = tempMicromorphicVarGrad;
    micromorphicStressGrad = tempMicromorphicStressGrad;    
}


void SecondGradientMaterialStatus :: initTempStatus()
//
// initialize record at the begining of new load step
//
{
    StructuralMaterialStatus :: initTempStatus();

    // see if vectors describing reached equilibrium are defined
    /* if ( this->giveStrainVector().giveSize() == 0 ) {
        strainVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
    }

    if ( this->giveStressVector().giveSize() == 0 ) {
        stressVector.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
	}*/

    // reset temp vars.
    tempMicromorphicVar = micromorphicVar;
    tempMicromorphicVarGrad = micromorphicVarGrad;
    tempMicromorphicStressGrad =  micromorphicStressGrad;
}


contextIOResultType
SecondGradientMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full ms context (saves state variables, that completely describe
// current state)
{
    contextIOResultType iores;

    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}


contextIOResultType
SecondGradientMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full material context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}

  /*
void MicromorphicMaterialStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const MicromorphicMaterialStatus &structStatus = dynamic_cast< MicromorphicMaterialStatus & >(tmpStat);

    strainVector = structStatus.giveStrainVector();
    stressVector = structStatus.giveStressVector();
    tempStressVector = structStatus.giveTempStressVector();
    tempStrainVector = structStatus.giveTempStrainVector();

    PVector = structStatus.givePVector();
    tempPVector = structStatus.giveTempPVector();
    CVector = structStatus.giveCVector();
    tempCVector = structStatus.giveTempCVector();
    FVector = structStatus.giveFVector();
    tempFVector = structStatus.giveTempFVector();
}

void MicromorphicMaterialStatus :: addStateVariables(const MaterialStatus &iStatus)
{
    printf("Entering MicromorphicMaterialStatus :: addStateVariables().\n");
}
  */
} // end namespace oofem

