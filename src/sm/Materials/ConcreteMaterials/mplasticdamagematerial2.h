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

#ifndef mplasticdamagematerial2_h
#define mplasticdamagematerial2_h

#include "../sm/Materials/structuralmaterial.h"
#include "Materials/linearelasticmaterial.h"
#include "dictionary.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "../sm/Materials/structuralms.h"
#include "Materials/ConcreteMaterials/mplasticmaterial2.h"

#include <vector>
#include <set>

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to MPlasticMaterial.
 * It is attribute of matStatusDictionary at every GaussPoint, for which this material
 * is active.
 *
 * Description:
 * Idea used there is that we have variables
 * describing:
 * -# state at previous equilibrium state (variables without temp)
 * -# state during searching new equilibrium (variables with temp)
 *    when we start search new state from previous equilibrium one we copy
 *    non-tem variables into temp ones. And after we reach new equilibrium
 *    (now described by temp variables) we copy tem-var into non-tepm ones
 *    (see function updateYourself).
 */
class MPlasticDamageMaterial2Status : public MPlasticMaterial2Status
{
public:

protected:
    /// Isotropic damage variables
    double damage, tempDamage;
public:
    MPlasticDamageMaterial2Status(int n, Domain * d, GaussPoint * g, int statusSize);
    virtual ~MPlasticDamageMaterial2Status();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    void letTempDamageBe(double v) { tempDamage = v; }
    double giveDamage() { return damage; }
    double giveTempDamage() { return tempDamage; }
    // definition
    virtual const char *giveClassName() const { return "MPlasticDamageMaterial2Status"; }
};

/**
 * This class represents a base class for non-associated multisurface plasticity.
 * The Multisurface plasticity is characterized by the following:
 * Let @f$\sigma, \varepsilon @f$, and @f$\varepsilon^p @f$ be the stress, total strain, and plastic strain vectors, respectively.
 * It is assumed that the total strain is decomposed into reversible elastic and irreversible plastic parts
 * @f$\varepsilon=\varepsilon^e+\varepsilon^p @f$.
 * The elastic response is characterized in terms of elastic constitutive matrix @f$ D^e @f$ as @f$ \sigma=D^e(\varepsilon-\varepsilon^e) @f$
 * As long as the stress remains inside the elastic domain, the deformation process is purely elastic and the
 * plastic strain does not change.
 *
 * It is assumed that the elastic domain, denoted as @f$ IE @f$ is bounded by a composite yield surface. It is defined as
 * @f[
 * IE=\{(\sigma,\kappa)|f_i(\sigma,\kappa)<0,\; i\in\{1,\cdots,m\}\}
 * @f]
 * where @f$ f_i(\sigma,\kappa) @f$ are @f$ m\ge1 @f$ yield functions intersecting in a possibly non-smooth fashion. The
 * vector @f$ \kappa @f$ contains internal variables controlling the evolution of yield surfaces (amount of hardening or softening).
 * The evolution of plastic strain @f$ \varepsilon^p @f$ is expressed in Koiter's form. Assuming the non-associated plasticity, this reads
 * @f[
 * \label{epe}
 * \varepsilon^p=\sum^{m}_{i=1} \lambda^i \partial_{\sigma}g_i(\sigma,\kappa)
 * @f]
 * where @f$ g_i @f$ are plastic potential functions. The @f$ \lambda^i @f$ are referred as plastic consistency parameters, which satisfy the following Kuhn-Tucker conditions
 * @f[
 * \label{ktc}
 * \lambda^i\ge0,\;f_i\le0,\;{\rm and}\ \lambda^i f_i=0
 * @f]
 * These conditions imply that in the elastic regime the yield function must remain negative and the rate of the plastic multiplier is zero
 * (plastic strain remains constant) while in the plastic regime the yield function must be equal to zero (stress remains on the surface) and the rate of the plastic multiplier is positive.
 * The evolution of vector of internal hardening/softening variables @f$ \kappa @f$  is expressed in terms of a general
 * hardening/softening law of the form
 * @f[
 * \dot{\kappa} = \dot{\kappa}(\sigma, \lambda)
 * @f]
 * where @f$ \lambda @f$ is the vector of plastic consistency parameters @f$ \lambda_i @f$.
 *
 */
class MPlasticDamageMaterial2 : public MPlasticMaterial2
{
protected:

public:
    MPlasticDamageMaterial2(int n, Domain * d);
    virtual ~MPlasticDamageMaterial2();

    virtual const char *giveClassName() const { return "MPlasticDamageMaterial2"; }

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }


    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);


    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *,
                                      const FloatArray &, TimeStep *);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }

    virtual double computeDamage(GaussPoint *gp, const FloatArray &strainSpaceHardeningVariables, TimeStep *tStep) = 0;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:


    void give_dLambda_dEps_Matrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    
    // next functions overloaded rom structural material level
    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    
};
} // end namespace oofem
#endif // mplasticdamagematerial2_h
