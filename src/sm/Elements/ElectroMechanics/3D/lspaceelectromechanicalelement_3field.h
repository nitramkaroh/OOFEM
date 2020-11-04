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

#ifndef lspaceelectromechanicalelement_h
#define lspaceelectromechanicalelement_h

#include "../sm/Elements/ElectroMechanics/baseelectromechanicalelement_3fields.h"
#include "Elements/3D/lspace.h"

#define _IFT_LSpaceElectroMechanicalElement_3fields_Name "lspaceelmechelem_3fields"

namespace oofem {
class FEI3dHexaLin;

/**
 * This class implements a Linear 3d 8-node finite element for stress analysis.
 * Each node has 3 degrees of freedom.
 *
 * One single additional attribute is needed for Gauss integration purpose :
 * 'jacobianMatrix'. This 3x3 matrix contains polynomials.
 * TASKS :
 * - Calculating its Gauss points.
 * - Calculating its B,D,N matrices and dV.
 */
class LSpaceElectroMechanicalElement  :  public LSpace, public BaseElectroMechanicalElement_3Fields
{
protected:
    static FEI3dHexaLin interpolation;
public:
    LSpaceElectroMechanicalElement_3Fields(int n, Domain * d);
    virtual ~LSpaceElectroMechanicalElement_3Fields() { }
    virtual FEInterpolation *giveInterpolation() const;

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LSpaceElectroMechanicalElement_3Fields_Name; }
    virtual const char *giveClassName() const { return "LSpaceElectroMechanicalElement_3Fields"; }
    virtual void postInitialize();

    virtual void computeElectricPotentialBmatrixAt(GaussPoint *gp, FloatMatrix &Be);
    virtual void computeElectricDisplacementNmatrixAt(GaussPoint *gp, FloatMatrix &Nd);

    virtual NLStructuralElement *giveStructuralElement(){return this;}

    virtual int giveNumberOfElectricPotentialDofs(){return 8;}
    virtual int giveNumberOfDisplacementDofs(){return 24;}
    virtual int giveNumberOfElectricDisplecementDofs(){return 24};
    virtual int giveNumberOfDofs()return{56;}


    void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer);
    virtual void giveDofManDofIDMask_phi(IntArray &answer);
    virtual void giveDofManDofIDMask_d(IntArray &answer);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep){BaseElectroMechanicalElement :: computeStiffnessMatrix(answer, mode, tStep);}
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord){BaseElectroMechanicalElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);}
    
    
};
} // end namespace oofem
#endif // lspaceelectromechanicalelement_h
