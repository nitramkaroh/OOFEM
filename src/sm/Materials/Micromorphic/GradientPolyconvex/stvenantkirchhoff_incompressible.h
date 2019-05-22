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

#ifndef stvenantkirchhoff_incompressible_h
#define stvenantkirchhoff_incompressible_h


#include "../sm/Materials/isolinearelasticmaterial.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "cltypes.h"


///@name Input fields for MicromorphLEmat
//@{
#define _IFT_StVenantKirchhoffIncompressibleMaterial_Name "stvenantkirchhoffincompressiblematerial"
#define _IFT_StVenantKirchhoffIncompressibleMaterial_K "k"
//@}

namespace oofem {


/**
 */
class StVenantKirchhoffIncompressibleMaterial : public IsotropicLinearElasticMaterial
{
protected:
  double K;
public:
    StVenantKirchhoffIncompressibleMaterial(int n, Domain * d);
    virtual ~StVenantKirchhoffIncompressibleMaterial();

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_StVenantKirchhoffIncompressibleMaterial_Name; }
    virtual const char *giveClassName() const { return "StVenantKirchhoffGradientPolyconvexMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);
    void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
  };
 

} // end namespace oofem
#endif
