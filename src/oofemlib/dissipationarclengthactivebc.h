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

#ifndef dissipationarclengthactivebc_h
#define dissipationarclengthactivebc_h

#include "arclengthactivebc.h"
#include "dofmanager.h"
#include <memory>

#define _IFT_DissipationArcLengthActiveBoundaryCondition_Name   "dissipationarclengthbc"

///@name Input fields for active boundary condition
//@{
#define _IFT_DissipationArcLengthActiveBoundaryCondition_dissipationLength   "deltal"
#define _IFT_DissipationArcLengthActiveBoundaryCondition_elementSet   "elemset"
//@}


namespace oofem {
/**
 * Abstract base class for arc lenght active boundary conditions.
 * This boundary condition specify the arc length constrain, i.e., it can have the classical form
 * (\Delta u)^T\Delta u  + \Delta \lambda \Psi \bar{f}_{ext}^t\bar{f}_{ext} = l_{arc}
 * or another form, like the dissipation based one, i.e., \Delta D_{diss} = l_{diss}
 */
class OOFEM_EXPORT DissipationArcLengthActiveBoundaryCondition : public ArcLengthActiveBoundaryCondition
{
 protected:
    int elementSet;
    double dl;
    std :: unique_ptr< DofManager > lm;
    IntArray refLoadLocationArray;
    IntArray elementList;
  
public:
    DissipationArcLengthActiveBoundaryCondition(int n, Domain * d);
    /// Destructor.
    virtual ~DissipationArcLengthActiveBoundaryCondition() { }

    double giveDeltaLambda() override {return 0;}

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms = NULL) override;
    virtual const char *giveClassName() const { return "DissipationArcLengthActiveBoundaryCondition"; }
    virtual const char *giveInputRecordName() const { return _IFT_DissipationArcLengthActiveBoundaryCondition_Name; }
    int giveNumberOfInternalDofManagers() override { return 1; }
    virtual DofManager *giveInternalDofManager(int i) { return this->lm.get(); }
    void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);
    void postInitialize() override;

 protected:
    void assemble(SparseMtrx &answer, TimeStep *tStep, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;
    void assembleReferenceLoadVector(SparseMtrx &answer, const IntArray &lambdaeq);


};
} // end namespace oofem
#endif // dissipationarclengthactivebc_h
