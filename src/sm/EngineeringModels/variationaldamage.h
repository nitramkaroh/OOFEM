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

#ifndef variationaldamage_h
#define variationaldamage_h

#include "../sm/EngineeringModels/nlinearstatic.h"
#include "sparsenonlinsystemnm.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"

///@name Input fields for NonLinearStatic
//@{
#define _IFT_VariationalDamage_Name "variationaldamage"
#define _IFT_VariationalDamage_MaxActivatedNodes "maxactivatednodes"
//@}

namespace oofem {

/**
 * This class implements nonlinear static engineering problem.
 **/
class VariationalDamage : public NonLinearStatic
{
 private:

    IntArray damageIndicatorArray;
    // Lists for each dof group
    std :: vector< CustomEquationNumbering > UnknownNumberingSchemeList;
    std :: vector< std :: unique_ptr< SparseMtrx > > stiffnessMatrixList;
    std :: vector< IntArray > locArrayList;
    void giveTotalLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s, Domain *d);
    int maxActivNodes;

protected:

public:
    VariationalDamage(int i, EngngModel * _master = NULL);
    virtual ~VariationalDamage();
    IRResultType initializeFrom(InputRecord *ir);
    void initializeYourself(TimeStep *tStep);
    //7jyu    NumericalMethod* giveNumericalMethod(MetaStep *mStep);
    void updateYourself(TimeStep *tStep);
    void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    virtual const char *giveClassName() const { return "StaggeredSolver"; }
    virtual const char *giveInputRecordName() const { return _IFT_VariationalDamage_Name; }
    void giveInternalForces(FloatArray &answer, bool normFlag, int di, TimeStep *tStep);
    
protected:
    void assemble(SparseMtrx &answer, TimeStep *tStep, const MatrixAssembler &ma,
		  const UnknownNumberingScheme &s, Domain *domain);

};
} // end namespace oofem
#endif // variationaldamage_h
