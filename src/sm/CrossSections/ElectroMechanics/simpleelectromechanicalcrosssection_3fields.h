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

#ifndef simpleelectromechanicalcrosssection_3fields_h
#define simpleelectromechanicalcrosssection_3fields_h

#include "crosssection.h"
#include "../sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for SimpleCrossSection
//@{
#define _IFT_SimpleElectroMechanicalCrossSection_3Fields_Name "simpleemcs"
#define _IFT_SimpleElectroMechanicalCrossSection_3Fields_thick "thick"
#define _IFT_SimpleElectroMechanicalCrossSection_3Fields_width "width"
#define _IFT_SimpleElectroMechanicalCrossSection_3Fields_area "area"
#define _IFT_SimpleElectroMechanicalCrossSection_3Fields_MaterialNumber "material" ///< Material number for the bulk material
//@}

namespace oofem {
/**
 * Class implementing "simple" electo-mechanical cross section model in finite element problem.
 * A cross section  is an attribute of a domain. It is usually also attribute of many
 * elements.
 *
 * The simple cross section implementation does not perform any integration over cross-section volume,
 * it represents a cross section model, where the whole cross section model is represented by single integration point.
 * and therefore all requests for characteristic contributions (stiffness) and for real stress computations are simply
 * passed to parent StructuralCrossSection class, which invokes corresponding material mode services.
 * Please note, that it is assumed that material model will support these material modes and provide
 * corresponding services for characteristic components and stress evaluation.
 * For description, how to incorporate more elaborate models of cross section, please read
 * base CrossSection documentation.
 *
 */
class OOFEM_EXPORT SimpleElectroMechanicalCrossSection_3Fields : public CrossSection
{
public:
    /**
     * Constructor.
     * @param n Cross section number.
     * @param d Associated domain.
     */
    SimpleElectroMechanicalCrossSection_3Fields(int n, Domain *d) : CrossSection(n, d) {
        materialNumber = 0;
        czMaterialNumber = 0;
    }


    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode mode);

    void give_FirstPKStressVector_ElectricFieldVector(FloatArray &P, FloatArray &E, GaussPoint *gp, const FloatArray &F, const FloatArray &D, TimeStep *tStep);
    void giveElectroMechanicalConstitutiveMatrix_dPdF(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp,TimeStep *tStep);
    void giveElectroMechanicalConstitutiveMatrix_dPdD(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp,TimeStep *tStep);
    void giveElectroMechanicalConstitutiveMatrix_dEdD(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp,TimeStep *tStep);
   
	
    

    /**
     * Initializes receiver acording to object description stored in input record.
     * Calls CrossSection initializeFrom service and reads the values of
     * - 'thick' thickness
     * - 'width' width
     * - 'area' area
     * @param ir Record to read off.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void createMaterialStatus(GaussPoint &iGP); // ES


    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "SimpleElectroMechanicalCrossSection_3Fields"; }
    virtual const char *giveInputRecordName() const { return _IFT_SimpleElectroMechanicalCrossSection_3Fields_Name; }

    virtual double give(int aProperty, GaussPoint *gp);
    virtual double give(CrossSectionProperty a, GaussPoint *gp) { return CrossSection :: give(a, gp); }
    virtual double give(CrossSectionProperty a, const FloatArray &coords, Element *elem, bool local) { return CrossSection :: give(a, coords, elem, local); }
    virtual int giveIPValue(FloatArray &answer, GaussPoint *ip, InternalStateType type, TimeStep *tStep);
    virtual Material *giveMaterial(IntegrationPoint *ip);

    int giveMaterialNumber() const { return this->materialNumber; }
    void setMaterialNumber(int matNum) { this->materialNumber = matNum; }
    virtual int checkConsistency();
    virtual Interface *giveMaterialInterface(InterfaceType t, IntegrationPoint *ip);



    virtual int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp);
    virtual int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp);
    virtual int estimatePackSize(DataStream &buff, GaussPoint *gp);

protected:
    int materialNumber;   // material number
    int czMaterialNumber; // cohesive zone material number
};
} // end namespace oofem
#endif // simpleelectromechanicalcrosssection_3fields_h
