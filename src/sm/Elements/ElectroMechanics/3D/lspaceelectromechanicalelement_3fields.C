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

#include "Elements/ElectroMechanics/3D/lspaceelectromechanicalelement_3fields.h"
#include "fei3dhexalin.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "crosssection.h"
#include "classfactory.h"



namespace oofem {
REGISTER_Element(LSpaceElectroMechanicalElement_3Fields);

FEI3dHexaLin LSpaceElectroMechanicalElement_3Fields :: interpolation;

  LSpaceElectroMechanicalElement_3Fields :: LSpaceElectroMechanicalElement_3Fields(int n, Domain *domain) : LSpace(n, domain), BaseElectroMechanicalElement_3Fields(n, domain)
    // Constructor.
{
}


FEInterpolation *LSpaceElectroMechanicalElement_3Fields :: giveInterpolation() const { return & interpolation; }



void
LSpaceElectroMechanicalElement_3Fields :: postInitialize()
{
    BaseElectroMechanicalElement_3Fields :: postInitialize();
    LSpace :: postInitialize();
    
}



void
LSpaceElectroMechanicalElement_3Fields :: computeElectricPotentialBmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3, this->giveNumberOfElectricPotentialDofs());
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1,  i ) = dNdx.at(i, 1);
        answer.at(2,  i ) = dNdx.at(i, 2);
	answer.at(3,  i ) = dNdx.at(i, 3);
    }

    answer.times(-1.);
}


void
LSpaceElectroMechanicalElement_3Fields :: computeElectricDisplacementNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatArray N; 
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(N, 3);
}
  

void
LSpaceElectroMechanicalElement_3Fields :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {D_u, D_v, D_w, E_phi, E_D1, E_D2, E_D3};
}
  

void
LSpaceElectroMechanicalElement_3Fields :: giveDofManDofIDMask_u(IntArray &answer)
{
  answer = {D_u, D_v, D_w};
}


void
LSpaceElectroMechanicalElement_3Fields :: giveDofManDofIDMask_phi(IntArray &answer)
{

  answer = {E_phi};

}

void
LSpaceElectroMechanicalElement_3Fields :: giveDofManDofIDMask_d(IntArray &answer)
{
  answer = {E_D1, E_D2, E_D3};
}

  

} // end namespace oofem
