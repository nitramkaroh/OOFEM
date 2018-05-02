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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

#include "../sm/Elements/Micromorphic/Microplastic/qplanestrainmicroplastic.h"
#include "fei2dquadquad.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(QPlaneStrainMicroplastic);

FEI2dQuadQuad QPlaneStrainMicroplastic :: interpolation(1, 2);

QPlaneStrainMicroplastic :: QPlaneStrainMicroplastic(int n, Domain *aDomain) : QPlaneStrain(n, aDomain), BaseMicromorphicElement()
    // Constructor.
{
  displacementDofsOrdering = {1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23};
  micromorphicDofsOrdering = {3,6,9,12,15,18,21,24};
}


void 
QPlaneStrainMicroplastic :: computeMicromorphicNMatrixAt(GaussPoint *gp,FloatMatrix &answer)
{
    FloatArray n;
    this->interpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 1);
}


void
QPlaneStrainMicroplastic :: computeMicromorphicBMatrixAt(GaussPoint *gp, FloatMatrix &answer)         
{
    FloatMatrix dNdx; 
    this->interpolation.evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beTranspositionOf(dNdx);
}


void
QPlaneStrainMicroplastic :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, M_MP};
}





void
QPlaneStrainMicroplastic :: giveDofManDofIDMask_u(IntArray &answer)
{
  answer = {D_u, D_v};
}


void
QPlaneStrainMicroplastic :: giveDofManDofIDMask_m(IntArray &answer)
{
  answer = {M_MP};
}

void
QPlaneStrainMicroplastic ::  postInitialize() 
{
  BaseMicromorphicElement :: postInitialize();
  QPlaneStrain :: postInitialize();
}

}
