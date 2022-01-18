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
REGISTER_Element(LSpaceElectroMechanicalElement_3Fields_D0);

FEI3dHexaLin LSpaceElectroMechanicalElement_3Fields_D0 :: interpolation;

  LSpaceElectroMechanicalElement_3Fields_D0 :: LSpaceElectroMechanicalElement_3Fields_D0(int n, Domain *domain) : LSpace(n, domain), BaseElectroMechanicalElement(n, domain)
    // Constructor.
{

  this->D0Node.reset( new ElementDofManager(1, aDomain, this) );
  this->D0Node->appendDof( new MasterDof(this->D0Node.get(), E_D1) );
  this->D0Node->appendDof( new MasterDof(this->D0Node.get(), E_D2) );
  this->D0Node->appendDof( new MasterDof(this->D0Node.get(), E_D3) );

}


FEInterpolation *LSpaceElectroMechanicalElement_3Fields_D0 :: giveInterpolation() const { return & interpolation; }



void
LSpaceElectroMechanicalElement_3Fields_D0 :: postInitialize()
{
    BaseElectroMechanicalElement_3Fields :: postInitialize();
    LSpace :: postInitialize();


    delete dofBCmap;
    dofBCmap = NULL;

    mBC.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, mBC, _IFT_DofManager_bc);
    bool hasBc = mBC.giveSize() != 0;
    IntArray dofIDArry(3);
    dofIDArry
    ///@todo This should eventually be removed, still here to preserve backwards compatibility:
    // check sizes
    if ( hasBc ) {
        if ( mBC.giveSize() != dofIDArry.giveSize() ) {
            OOFEM_ERROR("bc size mismatch. Size is %d and need %d", mBC.giveSize(), dofIDArry.giveSize());
        }
        this->dofBCmap = new std :: map< int, int >();
        for ( int i = 1; i <= mBC.giveSize(); ++i ) {
            if ( mBC.at(i) > 0 ) {
                ( * this->dofBCmap ) [ dofIDArry.at(i) ] = mBC.at(i);
            }
        }
    }
    
}



void
LSpaceElectroMechanicalElement_3Fields_D0 :: computeElectricFieldBmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3, this->giveNumberOfElectricDofs());
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1,  i ) = dNdx.at(i, 1);
        answer.at(2,  i ) = dNdx.at(i, 2);
	answer.at(3,  i ) = dNdx.at(i, 3);
    }

}


void
LSpaceElectroMechanicalElement_3Fields_D0 :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
  answer = {D_u, D_v, D_w, E_phi, E_D1, E_D2, E_D3};
}
  

void
LSpaceElectroMechanicalElement_3Fields_D0 :: giveDofManDofIDMask_u(IntArray &answer)
{
  answer = {D_u, D_v, D_w};
}


void
LSpaceElectroMechanicalElement_3Fields_D0 :: giveDofManDofIDMask_e(IntArray &answer)
{

  answer = {E_phi};

}

void
LSpaceElectroMechanicalElement_3Fields_D0 :: giveDofManDofIDMask_d(IntArray &answer)
{

  answer = {E_D1, E_D2, E_D3};

}

  

} // end namespace oofem
