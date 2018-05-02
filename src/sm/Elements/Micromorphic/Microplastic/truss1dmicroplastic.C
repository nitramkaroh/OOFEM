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

#include "../sm/Elements/Micromorphic/Microplastic/truss1dmicroplastic.h"
#include "fei1dlin.h"
#include "fei1dquad.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Truss1dMicroplastic);


Truss1dMicroplastic :: Truss1dMicroplastic(int n, Domain *aDomain) : Truss1d(n, aDomain), BaseMicromorphicElement()
    // Constructor.
{

  displacementDofsOrdering = {1,3};
  micromorphicDofsOrdering = {2,4};
  
}


void
Truss1dMicroplastic :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {D_u, M_MP};

}

void 
Truss1dMicroplastic :: computeMicromorphicNMatrixAt(GaussPoint *gp,FloatMatrix &answer)
{
  FloatArray n;
  this->interp.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
  answer.beNMatrixOf(n, 1);
}


void
Truss1dMicroplastic :: computeMicromorphicBMatrixAt(GaussPoint *gp, FloatMatrix &answer)         
{
    FloatMatrix dNdx; 
    interp.evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beTranspositionOf(dNdx);
}


void
Truss1dMicroplastic :: giveDofManDofIDMask_u(IntArray &answer)
{
  answer = {D_u};
}


void
Truss1dMicroplastic :: giveDofManDofIDMask_m(IntArray &answer)
{
  answer = {M_MP};
}


void
Truss1dMicroplastic ::  postInitialize() 
{
  BaseMicromorphicElement :: postInitialize();
  Truss1d :: postInitialize();
}

}
