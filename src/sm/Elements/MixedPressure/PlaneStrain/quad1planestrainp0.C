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

#include "Elements/MixedPressure/PlaneStrain/quad1planestrainp0.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"
#include "fei2dquadconst.h"
#include "masterdof.h"

namespace oofem {
REGISTER_Element(Quad1PlaneStrainP0);

FEI2dQuadConst Quad1PlaneStrainP0 :: interpolation(1, 2);

Quad1PlaneStrainP0 :: Quad1PlaneStrainP0(int n, Domain *aDomain) :Quad1PlaneStrain(n, aDomain)
{
  displacementDofsOrdering = {1,2,3,4,5,7,8};
  pressureDofsOrdering = {9};
  this->pressureNode.reset( new ElementDofManager(1, aDomain, this) );
  this->pressureNode->appendDof( new MasterDof(this->pressureNode.get(), P_f) );


}



 
void 
Quad1PlaneStrainP0 :: computePressureNMatrixAt(GaussPoint *gp,FloatArray &answer)
{
    FloatArray n;
    this->interpolation.evalN( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


void
Quad1PlaneStrainP0 :: computeVolumetricBmatrixAt(GaussPoint *gp, FloatArray &answer, NLStructuralElement *elem)
{
    answer.resize(8);
    FloatMatrix dN;
    elem->giveInterpolation()->evaldNdx( dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int j = 0, k = 0; j < 4; j++, k += 2 ) {
            answer(k)     = dN(j, 0);
            answer(k + 1) = dN(j, 1);
        }
}
   

void
Quad1PlaneStrainP0 :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {D_u, D_v};
}


void
Quad1PlaneStrainP0 :: giveDofManDofIDMask_u(IntArray &answer)
{
  answer = {D_u, D_v};
}


void
Quad1PlaneStrainP0 :: giveDofManDofIDMask_p(IntArray &answer)
{
  answer = {P_f};
}


void Quad1PlaneStrainP0 :: giveInternalDofManDofIDMask(int i, IntArray &answer) const
{
  answer = {P_f};
}

  

  
  
DofManager *
Quad1PlaneStrainP0 :: giveInternalDofManager(int i) const
{
  return pressureNode.get();
}


   
  
void
Quad1PlaneStrainP0 ::  postInitialize() 
{
  BaseMixedPressureElement :: postInitialize();
  Quad1PlaneStrain :: postInitialize();
}

  
  

} // end namespace oofem
