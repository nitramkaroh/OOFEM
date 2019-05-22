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

#ifndef hyperelasticmaterialextensioninterface_h
#define hyperelasticmaterialextensioninterface_h

#include "interface.h"
#include "domain.h"
#include "floatmatrix.h"
#include "floatarray.h"

///@name Input fields for HyperElasticMaterial
//@{
#define _IFT_HyperElasticMaterial_Name "hyperelmat"
#define _IFT_HyperElasticMaterial_k "k"
#define _IFT_HyperElasticMaterial_g "g"
//@}

namespace oofem {
/**
 * Saint Venant–Kirchhoff model defined by shear and bulk modulus.
 * @todo Should we even have this? Isn't this identical to the isotropic linear elastic model? / Mikael
 */
class HyperElasticMaterialExtensionInterface : public Interface
{
protected:
    
public:
    HyperElasticMaterialExtensionInterface(int n, Domain * d);

    /**
     * Method for computing the first invariant of right Cauchy-Green tensor C from deformation gradient F.
     * @return the First Invariant of C
     * @param F Deformation Gradient in the matrix form
     */
    virtual double compute_I1_C_from_F(const FloatMatrix &F);

    /**
     * Method for computing the second invariant of right Cauchy-Green tensor C from deformation gradient F.
     * @return the Second Invariant of C
     * @param F Deformation Gradient in the matrix form
     */
    virtual double compute_I2_C_from_F(const FloatMatrix &F);

    /**
     * Method for computing the third invariant of right Cauchy-Green tensor C from deformation gradient F.
     * @return the Third Invariant of C
     * @param F Deformation Gradient in the matrix form
     */
    virtual double compute_I3_C_from_F(const FloatMatrix &F);

    /**
     * Method for computing the first derivative of the first invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I1(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */
    virtual void compute_dI1dF(FloatArray &answer, const FloatMatrix &F);
    /**
     * Method for computing the first derivative of the second invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I2(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */
    virtual void compute_dI2dF(FloatArray &answer, const FloatMatrix &F);

     /**
     * Method for computing the first derivative of the third invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I3(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */
    virtual void compute_dI3dF(FloatArray &answer, const FloatMatrix &F);

    /**
     * Method for computing the second derivative of the first invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I1(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */
    virtual void compute_d2I1dF2(FloatMatrix &answer, const FloatMatrix &F);

    /**
     * Method for computing the second derivative of the second invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I1(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */

    virtual void compute_d2I2dF2(FloatMatrix &answer, const FloatMatrix &F);

    /**
     * Method for computing the second derivative of the third invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I1(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */
    virtual void compute_d2I3dF2(FloatMatrix &answer, const FloatMatrix &F);


    /**
     * Method for computing the first invariant of deviatoric part of right Cauchy-Green tensor C from deformation gradient F.
     * @param answer First Invariant of C
     * @param F Deformation Gradient in the matrix form
     */
    virtual double compute_I1_Cdev_from_F(const FloatMatrix &F);

    /**
     * Method for computing the second invariant of deviatoric part of right Cauchy-Green tensor C from deformation gradient F.
     * @param answer Second Invariant of C
     * @param F Deformation Gradient in the matrix form
     */
    virtual double compute_I2_Cdev_from_F(const FloatMatrix &F);

    /**
     * Method for computing the first derivative of the deviatoric part of the first invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I1(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */
    virtual void compute_dI1devdF(FloatArray &answer, const FloatMatrix &F);
    /**
     * Method for computing the first derivative of the deviatoric part of the second invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I2(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */
    virtual void compute_dI2devdF(FloatArray &answer, const FloatMatrix &F);

     
    /**
     * Method for computing the second derivative of the deviatoric part of the first invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I1(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */
    virtual void compute_d2I1devdF2(FloatMatrix &answer, const FloatMatrix &F);

    /**
     * Method for computing the second derivative of the deviatoric part of the second invariant of right Cauchy-Green tensor C wrt deformation gradient F.
     * @param answer Derivative of I1(C) wrt F
     * @param F Deformation Gradient in the matrix form
     */

    virtual void compute_d2I2devdF2(FloatMatrix &answer, const FloatMatrix &F);

    
    
    
};
} // end namespace oofem
#endif
