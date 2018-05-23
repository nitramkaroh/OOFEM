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

#ifndef pressurefollowerload_h
#define pressurefollowerload_h

#include "activebc.h"

#include <utility>
#include <list>

///@name Input fields for surface tension boundary condition
//@{
#define _IFT_PressureFollowerLoad_Name "pressurefollowerload"
#define _IFT_PressureFollowerLoad_pressure "pressure"
#define _IFT_PressureFollowerLoad_useTangent "usetangent"
//@}

namespace oofem {
class StructuralElement;
 class Element;

/**
 * Computes the load (and possibly tangent) for pressure follower load. 
 * Supports 3d elements, 2d memebrane element, and axisymetric 1d membrane element
 */
class OOFEM_EXPORT PressureFollowerLoad : public ActiveBoundaryCondition
{
    double pressure; ///< pressure.
    bool useTangent; ///< Determines if tangent should be used.

public:
    /**
     * Constructor. Creates boundary an active condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    PressureFollowerLoad(int n, Domain * d) : ActiveBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~PressureFollowerLoad() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void assemble(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorms = NULL);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual const char *giveClassName() const { return "PressureFollowerLoad"; }
    virtual const char *giveInputRecordName() const { return _IFT_PressureFollowerLoad_Name; }
    virtual int giveApproxOrder() { return 0; }
    void computeVolumeLoadVectorFromElement(FloatArray &answer, Element *e, int iSurf, TimeStep *tStep){;}
    

protected:
    /**
     * Helper function for computing the contributions to the load vector.
     */
    void computeLoadVectorFromElement(FloatArray &answer, Element *e, int side, TimeStep *tStep);
    /**
     * Helper function for computing the tangent (@f$ K = \frac{\mathrm{d}F}{\mathrm{d}u} @f$)
     */
    void computeTangentFromElement(FloatMatrix &answer, Element *e, int side, TimeStep *tStep);

    void giveSurfacedNdxi_dDdDMatrces(FloatMatrix &dNk_dD, FloatMatrix &dNe_dD, const FloatMatrix &dNdxi, int nNodes);
};
} // end namespace oofem
#endif // pressurefollowerload_h
