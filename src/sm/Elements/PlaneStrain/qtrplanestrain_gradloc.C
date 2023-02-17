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

#include "../sm/Elements/PlaneStrain/qtrplanestrain_gradloc.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "classfactory.h"
#include "fei2dtrquad.h"


namespace oofem {
REGISTER_Element(QTrPlaneStrainGradLoc);


QTrPlaneStrainGradLoc :: QTrPlaneStrainGradLoc(int n, Domain *aDomain) :
  QTrPlaneStrain(n, aDomain), BaseSecondGradientElement()
{

}



void
QTrPlaneStrainGradLoc :: computeGmatrixAt(GaussPoint *gp, FloatMatrix &answer)                
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix d2Ndx2; 
    interp->evald2Ndx2( d2Ndx2, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(8, this->giveElement()->computeNumberOfDofs());
    answer.zero();

    for ( int i = 1; i <= d2Ndx2.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = d2Ndx2.at(i, 1);
        answer.at(2, i * 2 - 1) = d2Ndx2.at(i, 3);

     	answer.at(3, i * 2 - 0) = d2Ndx2.at(i, 3);
        answer.at(4, i * 2 - 0) = d2Ndx2.at(i, 2);

        answer.at(5, i * 2 - 1) = d2Ndx2.at(i, 3);
        answer.at(6, i * 2 - 1) = d2Ndx2.at(i, 2);

        answer.at(7, i * 2 - 0) = d2Ndx2.at(i, 1);
        answer.at(8, i * 2 - 0) = d2Ndx2.at(i, 3);
    }
    
}
  
} // end namespace oofem
