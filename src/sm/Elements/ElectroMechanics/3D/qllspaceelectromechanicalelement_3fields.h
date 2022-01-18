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

#ifndef qllspaceelectromechanicalelement_3fields_h
#define qllspaceelectromechanicalelement_3fields_h

#include "../sm/Elements/ElectroMechanics/baseelectromechanicalelement_3fields.h"
#include "Elements/3D/qspace.h"

#define _IFT_QLLSpaceElectroMechanicalElement_3Fields_Name "qllspaceelmechelem_3fields"

namespace oofem {
class FEI3dHexaLin;
 class FEI3dHexaQuad;

/**
 * This class implements a Quadratic in u and phi and linead in D 3d 20-node finite element for stress analysis.
 *
 * One single additional attribute is needed for Gauss integration purpose :
 * 'jacobianMatrix'. This 3x3 matrix contains polynomials.
 * TASKS :
 * - Calculating its Gauss points.
 * - Calculating its B,D,N matrices and dV.
 */
class QLLSpaceElectroMechanicalElement_3Fields  :  public QSpace, public BaseElectroMechanicalElement_3Fields
{
protected:
    static FEI3dHexaLin interpolation_lin;
    static FEI3dHexaQuad interpolation;
public:
    QLLSpaceElectroMechanicalElement_3Fields(int n, Domain * d);
    virtual ~QLLSpaceElectroMechanicalElement_3Fields() { }
    virtual FEInterpolation *giveInterpolation() const;

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QLLSpaceElectroMechanicalElement_3Fields_Name; }
    virtual const char *giveClassName() const { return "QLLSpaceElectroMechanicalElement_3Fields"; }

    FEInterpolation* giveInterpolation_lin() const;

    virtual void computeElectricPotentialBmatrixAt(GaussPoint *gp, FloatMatrix &Be);
    virtual void computeElectricDisplacementNmatrixAt(GaussPoint *gp, FloatMatrix &Nd);

    virtual NLStructuralElement *giveStructuralElement(){return this;}

    virtual int giveNumberOfElectricPotentialDofs(){return 8;}
    virtual int giveNumberOfDisplacementDofs(){return 60;}


    
    /*    virtual int giveNumberOfElectricDisplacementDofs(){return 60;}
    virtual int giveNumberOfDofs() {return 140;}
    */
    virtual int giveNumberOfElectricDisplacementDofs(){return 24;}
    virtual int giveNumberOfDofs() {return 92;}


    void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer);
    virtual void giveDofManDofIDMask_phi(IntArray &answer);
    virtual void giveDofManDofIDMask_d(IntArray &answer);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep){BaseElectroMechanicalElement_3Fields :: computeStiffnessMatrix(answer, mode, tStep);}
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord){BaseElectroMechanicalElement_3Fields :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);}
 protected:
    void postInitialize() override;
    
};
} // end namespace oofem
#endif // lspaceelectromechanicalelement_h
