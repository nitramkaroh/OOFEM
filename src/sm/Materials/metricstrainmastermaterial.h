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

#ifndef metricstrainmastermaterial_h
#define metricstrainmastermaterial_h

#include "structuralmaterial.h"
#include "structuralms.h"
#include "linearelasticmaterial.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for LargeStrainMasterMaterial
//@{
#define _IFT_MetricStrainMasterMaterial_Name "msmastermat"
#define _IFT_MetricStrainMasterMaterial_m "m"
#define _IFT_MetricStrainMasterMaterial_slaveMat "slavemat"
//@}

namespace oofem {
class GaussPoint;
class Domain;

/**
 * Metric strain master material.
 * Stress wokr-conjugated to metric strain and corresponding stiffness are computed from small strain(slaveMat) material model using a strain tensor from the Metric strain family (depending on parameter m)
 * Choice m = 0 leads to a second order approximation of the logarithmic strain, m = 2 corresponds to the Green-Lagrangian strain ...)
 * Stress and Stiffness are converted to 1.PK stress and first elasticity stiffness
 * @Reference:
 * A. Curnier, Ph. Zysset
 * A family of metric strain and conjugate stresses, prolonging usual material laws from small to large              * transformations, International Journal of Solid and Structures, 43 (2006), 3057-3086
 * @author Martin Horák 
 * @date 2014-25-04
 * @todo Implementation of PlaneStress, PlaneStrain and 1D cases
 */
class MetricStrainMasterMaterial : public StructuralMaterial
{
protected:
    /// Reference to the basic elastic material.
    LinearElasticMaterial *linearElasticMaterial;

    /// 'slave' material model number.
    int slaveMat;
    /// Specifies the strain tensor.
    double m;


public:
    MetricStrainMasterMaterial(int n, Domain * d);
    virtual ~MetricStrainMasterMaterial();

    //    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual const char *giveInputRecordName() const { return _IFT_MetricStrainMasterMaterial_Name; }
    virtual const char *giveClassName() const { return "MetricStrainMasterMaterial"; }

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode,
                                               GaussPoint *gp,
                                               TimeStep *tStep){;}

    virtual void give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                               MatResponseMode,
                                               GaussPoint *gp,
						    TimeStep *tStep);



    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *, const FloatArray &, TimeStep *)
    { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    virtual void giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep);

    /// transformation matrices
    void  compute_Sm_2_P_TransformationMatrix(FloatMatrix &answer, FloatMatrix &F);
    void  compute_Sm_2_P_TransformationMatrix(FloatMatrix &answer, FloatMatrix &F, FloatMatrix &invF);
     void convert_dSmdEm_2_dPdF(FloatMatrix &answer, const FloatMatrix &dSdE, FloatArray &vSm, FloatArray &vF, MaterialMode matMode);  


    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
};

//=============================================================================


class MetricStrainMasterMaterialStatus : public StructuralMaterialStatus
{
protected:
  //    FloatMatrix Pmatrix, TLmatrix, transformationMatrix;
    int slaveMat;
    /// Gauss Point for slave material status
    GaussPoint *slaveGp;
    /// Slave material status
    StructuralMaterialStatus *slaveStatus;

public:
    MetricStrainMasterMaterialStatus(int n, Domain * d, GaussPoint * g, int s);
    virtual ~MetricStrainMasterMaterialStatus();   

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    void printMasterMaterialOutputAt(FILE *file, TimeStep *tStep);
    void printSlaveMaterialOutputAt(FILE *File, TimeStep *tNow);

    GaussPoint* giveSlaveGaussPoint()
    { return this->slaveGp;}
    MaterialStatus *giveSlaveStatus() const
    {return this->slaveStatus;}
    void createSlaveStatus();


    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "LargeStrainMasterMaterialStatus"; }
};
} // end namespace oofem
#endif // metricstrainmastermaterial_h
