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

#ifndef pressurefollowerloadinterface_h
#define pressurefollowerloadinterface_h

#include "interface.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "element.h"
#include "feinterpol.h"
#include "gausspoint.h"


namespace oofem {
class GaussPoint;
class IntegrationRule;
class TimeStep;
/**
 * Provides Pressure Follower Load for an element.
 * @author Martin Horak
 */
class OOFEM_EXPORT PressureFollowerLoadElementInterface : public Interface
{
public:
    Element *element;


    /// Constructor.
    PressureFollowerLoadElementInterface(Element *e);
    virtual ~PressureFollowerLoadElementInterface();

    virtual const char *giveClassName() const { return "PressureFollowerLoadInterface"; }
    //    std :: string errorInfo(const char *func) const { return std :: string( giveClassName() ) + func; }


    
    // private:
    virtual void  surfaceEvalNumericalStiffMatrixAt(FloatMatrix &answer, FloatArray &dxdeta, FloatArray &dxdxi, int iSurf, GaussPoint *gp, TimeStep *tStep){answer.zero();}
    virtual void surfaceEvalNmatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp) = 0;
    virtual void surfaceEvaldNdxi(FloatMatrix &answer, int iSurf, GaussPoint *gp) = 0;
    virtual void surfaceEvalDeformedNormalAt(FloatArray &answer, FloatArray &dxdksi, FloatArray &dxdeta, int iSurf, GaussPoint *gp, TimeStep *tStep) = 0;

    virtual void surfaceEvalNormalDerivative(FloatMatrix &answer, int iSurf, GaussPoint *gp, TimeStep *tStep){;}
    
    virtual IntegrationRule* surfaceGiveIntegrationRule(int order, int iSurf) = 0;

};
} // end namespace oofem
#endif // pressurefollowerloadinterface
