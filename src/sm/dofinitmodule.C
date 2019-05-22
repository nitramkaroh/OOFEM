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

// Initialization module reading data related to Gauss points from a specified file

#include "gpinitmodule.h"
#include "floatarray.h"
#include "domain.h"
#include "engngm.h"
#include "classfactory.h"
#include <cassert>

namespace oofem {
REGISTER_InitModule(DofManagerInitModule)

DofManagerInitModule :: DofManagerInitModule(int n, EngngModel *e) : InitModule(n, e)
{ }


DofManagerInitModule :: ~DofManagerInitModule()
{ }


IRResultType
DofManagerInitModule :: initializeFrom(InputRecord *ir)
{
    return InitModule :: initializeFrom(ir);

    values.clear();
    IR_GIVE_FIELD(ir, values, _IFT_DofManagerInitModule_values);

    dofs.clear();
    IR_GIVE_FIELD(ir, dofs, _IFT_DofManagerInitModule_dofs);
    this->dofidmask = new IntArray(dofs);

    set = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, set, _IFT_DofManagerInitModule_set);
    
}


void
DofManagerInitModule :: giveValue(FloatArray &answer, FloatArray &dofIdMask, int set)
{
  answer = this->values;
  dofIdMask = this->dofidmask;
  set = this->set;  
}
  
} // namespace oofem
