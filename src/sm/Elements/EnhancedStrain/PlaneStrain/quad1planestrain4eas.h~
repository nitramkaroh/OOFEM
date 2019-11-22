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
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef quad1planestrain4eas_h
#define quad1planestrain4eas_h

#include "quad1planestrain.h"
#include "enhancedassumestrainelementinterface.h"


#define _IFT_Quad1PlaneStrain4EAS_Name "quad1planestrain4eas"

namespace oofem {
/**
 * This class implements an isoparametric four-node quadrilateral plane-
 * strain structural finite element enhanced by 4 modes
 */
class Quad1PlaneStrain4EAS : public Quad1PlaneStrain, public EnhancedAssumedStrainElementExtensionInterface
{
protected:
    
public:
    Quad1PlaneStrain4EAS(int n, Domain *d);
	virtual ~Quad1PlaneStrain4EAS(){;}

    
	virtual Interface *giveInterface(InterfaceType t) { if ( t == EnhancedAssumedStrainElementExtensionInterfaceType ) { return static_cast< EnhancedAssumedStrainElementExtensionInterface * >( this ); } else { return NULL; } }
	        
	virtual void computeEnhancedBmatrixAt(GaussPoint *gp, FloatMatrix &answer, StructuralElement *elem);
	virtual void computeEnhancedBHmatrixAt(GaussPoint *gp, FloatMatrix &answer, NLStructuralElement *elem);
    virtual int computeNumberOfEnhancedDofs(){return 4;}

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Quad1PlaneStrain4EAS_Name; }
    virtual const char *giveClassName() const { return "Quad1PlaneStrain4EAS"; }
    virtual classType giveClassID() const { return Quad1PlaneStrain4EASClass; }
    //virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _PlaneStrain; }

protected:
	void computeEnhancedModesAt(std::vector<FloatMatrix> &answer, GaussPoint* gp);
};
} // end namespace oofem
#endif // quad1planestrain4eas_h
