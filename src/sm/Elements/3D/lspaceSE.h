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

#ifndef lspacese_h
#define lspacese_h

#include "Elements/3D/lspace.h"


#define _IFT_LSpaceSE_Name "lspacese"

namespace oofem {
class FEI3dHexaLin;

/**
Description
 */
class LSpaceSE  : public LSpace
{
protected:
  static FEI3dHexaLin interpolation;

public:
    LSpaceSE(int n, Domain * d);
    virtual ~LSpaceSE() { }

    void computeStiffnessMatrix(FloatMatrix &answer,MatResponseMode rMode, TimeStep *tStep);
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
 protected:

        virtual void computeNlBmatrixAt(GaussPoint *gp, FloatMatrix &answer,FloatMatrix &G, TimeStep *tStep = NULL, int = 0, int = ALL_STRAINS);

};
} // end namespace oofem
#endif // lspacese_h
