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

#ifndef qbrickgradpolyconvex_h
#define qbrickgradpolyconvex_h

#include "../sm/Elements/3D/qspace.h"
#include "../sm/Elements/Micromorphic/basemicromorphicelement.h"

#define _IFT_QLBrickGradPolyconvex_Name "qlbrickgradpolyconvex"

namespace oofem {
class FEI3dhexalin;

class QLBrickGradPolyconvex : public QSpace, public BaseMicromorphicElement
{
protected:
    static FEI3dHexaLin micro_interpolation;

public:
    QLBrickGradPolyconvex(int n, Domain * d);
    virtual ~QLBrickGradPolyconvex() { }

 protected:


    virtual void computeMicromorphicNMatrixAt(GaussPoint *gp, FloatMatrix &N);
    virtual void computeMicromorphicBMatrixAt(GaussPoint *gp, FloatMatrix &B);

    virtual NLStructuralElement *giveElement() { return this; }
 public:
    virtual const char *giveInputRecordName() const { return _IFT_QLBrickGradPolyconvex_Name; }
    virtual const char *giveClassName() const { return "LBrickGradPolyconvex"; }

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer);
    virtual void giveDofManDofIDMask_m(IntArray &answer);

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep){BaseMicromorphicElement :: computeStiffnessMatrix(answer, mode, tStep);}
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord){BaseMicromorphicElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);}
    

    virtual int giveNumberOfMicromorphicDofs(){return 72;}
    virtual int giveNumberOfDisplacementDofs(){return 60;}
    virtual int giveNumberOfDofs(){return 132;}
    virtual int computeNumberOfDofs() {return this->giveNumberOfDofs();}

    virtual void postInitialize();
 protected:
    virtual void computeMicromorphicVars(FloatArray &micromorphVar, FloatArray &micromorphVarGrad, IntArray IdMask_m, GaussPoint *gp, TimeStep *tStep);   

};
} // end namespace oofem
#endif // lbrickgradpolyconvex_h
