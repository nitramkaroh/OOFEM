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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef martensitematerial_h
#define martensitematerial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"

///@name Input fields for OgdenMaterial
//@{
#define _IFT_MartensiteMaterial_Name "martensitematerial"
#define _IFT_MartensiteMaterial_m "m"
#define _IFT_MartensiteMaterial_D11 "d11"
#define _IFT_MartensiteMaterial_D12 "d12"
#define _IFT_MartensiteMaterial_D44 "d44"
#define _IFT_MartensiteMaterial_eta1 "eta1"
#define _IFT_MartensiteMaterial_eta2 "eta2"
#define _IFT_MartensiteMaterial_M "m"
//@}

namespace oofem {
/**
 * This class implements Multi-well Martenzite material
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 */
class MartensiteMaterial : public StructuralMaterial
{
protected:
    // Material parameters

    double M, D11, D12, D44;
    FloatArray arrayEpsilon;
    std::vector<FloatMatrix> vector_iF;
    std::vector<FloatArray>  vector_Ci;
    


public:
    MartensiteMaterial(int n, Domain *d);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep);
    

    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep);


    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);

       virtual const char *giveInputRecordName() const { return _IFT_MartensiteMaterial_Name; }
    virtual const char *giveClassName() const { return "MartensiteMaterial"; }
 protected:
     
    void giveSecondPKStressVector_3d(FloatMatrix &answer, GaussPoint *gp, const FloatMatrix &C, TimeStep *tStep);
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);


};



/**
 * This class implements associated Material Status to IsotropicDamageMaterial1.
 * Stores the characteristic length of the element.
 */
class MartensiteMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Scalar measure of dissipation
    double dissipation;
    /// Non-equilibrated scalar measure of dissipation.
    double tempDissipation;
    /// Index of activated energy
    int index = 1;
    /// Dissipation vector
    FloatArray dissipationVector;
    /// Non-equilibrated dissipation vector.
    FloatArray tempDissipationVector;
    // array of distances form well
    FloatArray array_wellDistance;

    bool init = false;
  
public:
    /// Constructor
    MartensiteMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~MartensiteMaterialStatus() { }

    // definition
    virtual const char *giveClassName() const { return "MartensiteMaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);


    /// Returns equilibrated dissipation
    double giveDissipation() { return dissipation; }
    /// Sets nonequilibrated dissipation
    void setTempDissipation(double td) { tempDissipation = td; }
    double giveTempDissipation() { return tempDissipation; }


    /// Returns equilibrated dissipation
    FloatArray &giveDissipationVector() { return dissipationVector; }
    /// Sets nonequilibrated dissipation
    void setTempDissipationVector(const FloatArray &tdv) { tempDissipationVector = tdv; }
    void setDissipationVector(const FloatArray &dv) { dissipationVector = dv; }


    double giveActiveWellIndex(){return index;}
    void setActiveWellIndex(double i){index = i;}
    
    void setArrayWellDistance(const FloatArray &awd){array_wellDistance = awd;}
    FloatArray &giveArrayWellDistance(){return array_wellDistance;}
    bool isDissipationVectorInitialized(){return init;}
    void setDissipationVectorInitialized(){init = true;}
    
    
};




 
} // end namespace oofem
#endif
