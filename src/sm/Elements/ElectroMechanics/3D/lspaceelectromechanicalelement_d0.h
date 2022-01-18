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

#ifndef lspaceelectromechanicalelement_d0_h
#define lspaceelectromechanicalelement_d0_h

#include "../sm/Elements/ElectroMechanics/baseelectromechanicalelement.h"
#include "Elements/3D/lspace.h"
#include "../sm/CrossSections/ElectroMechanics/simpleelectromechanicalcrosssection_3fields.h"


#define _IFT_LSpaceElectroMechanicalElement_D0_Name "lspaceelmechelem_d0"
#define _IFT_LSpaceElectroMechanicalElement_D0_T "t"

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
class LSpaceElectroMechanicalElement_D0  :  public LSpace, public BaseElectroMechanicalElement
{
protected:
    static FEI3dHexaLin interpolation;
    FloatArray D;
    FloatArray deltaU;
    FloatArray deltaPhi;
    int t = 0;
    ///
    FloatArray fd;
    FloatMatrix ikdd, kud, kphid;
    void setInverseOfKdd(const FloatMatrix &k){ikdd = k;}
    void setKud(const FloatMatrix &k) {this->kud = k;}
    void setKphid(const FloatMatrix &k) {this->kphid = k;}
    void setFd(const FloatArray &f){this->fd = f;}


    void giveInverseOfKdd(FloatMatrix &answer){answer = this->ikdd;}
    void giveKud(FloatMatrix &answer) {answer = this->kud;}
    void giveKphid(FloatMatrix &answer) {answer = this->kphid;}
    void giveFd(FloatArray &answer){answer = this->fd;}

    
public:
    LSpaceElectroMechanicalElement_D0(int n, Domain * d);
    virtual ~LSpaceElectroMechanicalElement_D0() { }
    virtual FEInterpolation *giveInterpolation() const;
    //
    void updateYourself(TimeStep *tStep);
    IRResultType initializeFrom(InputRecord *ir);


    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_LSpaceElectroMechanicalElement_D0_Name; }
    virtual const char *giveClassName() const { return "LSpaceElectroMechanicalElement_D0"; }
    virtual void postInitialize();

    virtual void computeElectricFieldBmatrixAt(GaussPoint *gp, FloatMatrix &Be);
    void computeElectricDisplacementNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual NLStructuralElement *giveStructuralElement(){return this;}
    virtual int giveNumberOfElectricDofs(){return 8;}
    virtual int giveNumberOfDisplacementDofs(){return 24;}
    virtual int giveNumberOfDofs(){return 32;}
    SimpleElectroMechanicalCrossSection_3Fields *giveCrossSection(); 

    void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer);
    virtual void giveDofManDofIDMask_e(IntArray &answer);
    
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep);

    void computeStiffnessMatrices(FloatMatrix &answer_uu, FloatMatrix &answer_up, FloatMatrix &answer_pp, MatResponseMode rMode, TimeStep *tStep);

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    
    void computeElectricDisplacementVector(FloatArray &answer,GaussPoint *gp, TimeStep *tStep);
    void updateElectricDisplacementVector(TimeStep *tStep);

    void compute_FirstPKStressVector_ElectricFieldVector(FloatArray &P, FloatArray &E, GaussPoint *gp, TimeStep *tStep);
    void computeElectricPotentialGradientVector(FloatArray &answer,GaussPoint *gp, TimeStep *tStep);
    
};
} // end namespace oofem
#endif // lspaceelectromechanicalelement_h
