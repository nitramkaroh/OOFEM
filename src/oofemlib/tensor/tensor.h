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

#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "tensor/FTensor.hpp"

#pragma once

using namespace FTensor;

static FTensor::Index<'i', 3> i_3;
static FTensor::Index<'j', 3> j_3;
static FTensor::Index<'k', 3> k_3;
static FTensor::Index<'l', 3> l_3;
static FTensor::Index<'m', 3> m_3;
static FTensor::Index<'n', 3> n_3;
/*static FTensor::Index<'o', 3> o_3;
static FTensor::Index<'p', 3> p_3;
static FTensor::Index<'q', 3> q_3;
static FTensor::Index<'r', 3> r_3;
static FTensor::Index<'s', 3> s_3;
static FTensor::Index<'t', 3> t_3;
static FTensor::Index<'u', 3> u_3;
static FTensor::Index<'v', 3> v_3;
static FTensor::Index<'w', 3> w_3;
static FTensor::Index<'x', 3> x_3;
static FTensor::Index<'y', 3> y_3;
static FTensor::Index<'z', 3> z_3;
*/
namespace oofem {

  inline FloatArrayF<9> to_voigt_form(const Tensor2_3d &t)
  {
    return {
      t(0, 0),
	t(1, 1),
	t(2, 2),
	t(1, 2),
	t(0, 2),
	t(0, 1),
	t(2, 1), 
	t(2, 0),
	t(1, 0)
	};
  }


  inline FloatMatrixF<9,9> to_voigt_form(const Tensor4_3d &t)
    {
    return {t(0,0,0,0),t(1,1,0,0),t(2,2,0,0),t(1,2,0,0),t(0,2,0,0),t(0,1,0,0),t(2,1,0,0),t(2,0,0,0),t(1,0,0,0),t(0,0,1,1),t(1,1,1,1),t(2,2,1,1),t(1,2,1,1),t(0,2,1,1),t(0,1,1,1),t(2,1,1,1),t(2,0,1,1),t(1,0,1,1),t(0,0,2,2),t(1,1,2,2),t(2,2,2,2),t(1,2,2,2),t(0,2,2,2),t(0,1,2,2),t(2,1,2,2),t(2,0,2,2),t(1,0,2,2),t(0,0,1,2),t(1,1,1,2),t(2,2,1,2),t(1,2,1,2),t(0,2,1,2),t(0,1,1,2),t(2,1,1,2),t(2,0,1,2),t(1,0,1,2),t(0,0,0,2),t(1,1,0,2),t(2,2,0,2),t(1,2,0,2),t(0,2,0,2),t(0,1,0,2),t(2,1,0,2),t(2,0,0,2),t(1,0,0,2),t(0,0,0,1),t(1,1,0,1),t(2,2,0,1),t(1,2,0,1),t(0,2,0,1),t(0,1,0,1),t(2,1,0,1),t(2,0,0,1),t(1,0,0,1),t(0,0,2,1),t(1,1,2,1),t(2,2,2,1),t(1,2,2,1),t(0,2,2,1),t(0,1,2,1),t(2,1,2,1),t(2,0,2,1),t(1,0,2,1),t(0,0,2,0),t(1,1,2,0),t(2,2,2,0),t(1,2,2,0),t(0,2,2,0),t(0,1,2,0),t(2,1,2,0),t(2,0,2,0),t(1,0,2,0),t(0,0,1,0),t(1,1,1,0),t(2,2,1,0),t(1,2,1,0),t(0,2,1,0),t(0,1,1,0),t(2,1,1,0),t(2,0,1,0),t(1,0,1,0)};
  }



  inline Tensor2_3d compute_tensor_cross_product(const Tensor2_3d &A, const Tensor2_3d &B) 
  {
    Tensor2_3d C (A(1,1) * B(2,2) - A(1,2) * B(2,1) - A(2,1) * B(1,2) + A(2,2) * B(1,1), A(1,2)*B(2,0) - A(1,0)*B(2,2) + A(2,0)*B(1,2) - A(2,2)*B(1,0), A(1,0)*B(2,1) - A(1,1)*B(2,0) - A(2,0)*B(1,1) + A(2,1)*B(1,0), A(0,2)*B(2,1) - A(0,1)*B(2,2) + A(2,1)*B(0,2) - A(2,2)*B(0,1), A(0,0)*B(2,2) - A(0,2)*B(2,0) - A(2,0)*B(0,2) + A(2,2)*B(0,0),A(0,1)*B(2,0) - A(0,0)*B(2,1) + A(2,0)*B(0,1) - A(2,1)*B(0,0), A(0,1)*B(1,2) - A(0,2)*B(1,1) - A(1,1)*B(0,2) + A(1,2)*B(0,1), A(0,2)*B(1,0) - A(0,0)*B(1,2) + A(1,0)*B(0,2) - A(1,2)*B(0,0), A(0,0)*B(1,1) - A(0,1)*B(1,0) - A(1,0)*B(0,1) + A(1,1)*B(0,0));
    return C;    
  }


  inline Tensor4_3d compute_tensor_cross_product(const Tensor2_3d &A)  
  {
    Tensor4_3d Ax(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    Ax(1,1,0,0) =  A(2,2);
    Ax(1,2,0,0) = -A(2,1);
    Ax(2,1,0,0) = -A(1,2);
    Ax(2,2,0,0) =  A(1,1);
    
    Ax(0,1,1,0) = -A(2,2);
    Ax(0,2,1,0) =  A(2,1); 
    Ax(2,1,1,0) =  A(0,2);
    Ax(2,2,1,0) = -A(0,1);

    Ax(0,1,2,0) =  A(1,2);
    Ax(0,2,2,0) = -A(1,1);
    Ax(1,1,2,0) = -A(0,2);
    Ax(1,2,2,0) =  A(0,1);

    Ax(1,0,0,1) = -A(2,2);
    Ax(1,2,0,1) =  A(2,0);
    Ax(2,0,0,1) =  A(1,2); 
    Ax(2,2,0,1) = -A(1,0);
 
    Ax(0,0,1,1) =  A(2,2);
    Ax(0,2,1,1) = -A(2,0);
    Ax(2,0,1,1) = -A(0,2);
    Ax(2,2,1,1) =  A(0,0);
  
    Ax(0,0,2,1) = -A(1,2);
    Ax(0,2,2,1) =  A(1,0);
    Ax(1,0,2,1) =  A(0,2);
    Ax(1,2,2,1) = -A(0,0);

    Ax(1,0,0,2) =  A(2,1);
    Ax(1,1,0,2) = -A(2,0);
    Ax(2,0,0,2) = -A(1,1);
    Ax(2,1,0,2) =  A(1,0);

    Ax(0,0,1,2) = -A(2,1);   
    Ax(0,1,1,2) =  A(2,0);
    Ax(2,0,1,2) =  A(0,1);
    Ax(2,1,1,2) = -A(0,0);

    Ax(0,0,2,2) =  A(1,1);
    Ax(0,1,2,2) = -A(1,0);
    Ax(1,0,2,2) = -A(0,1);
    Ax(1,1,2,2) =  A(0,0);
  
    return Ax;    
  }
 
  
  inline Tensor2_3d compute_cofactor(const Tensor2_3d &F)
  {

    Tensor2_3d cofF;
    cofF(i_3, j_3) = 0.5 * compute_tensor_cross_product(F,F)(i_3, j_3);
    return cofF;
  }


  inline double compute_determinant(const Tensor2_3d &F)
  {
    return ( 1./6. *  compute_tensor_cross_product(F,F)(m_3,n_3) * F(m_3,n_3) );
  }

 


  inline std::pair<double, Tensor2_3d>  compute_determinant_and_cofactor(const Tensor2_3d &F) 
    {

      Tensor2_3d cofF;
      cofF(i_3,j_3)=  0.5 * compute_tensor_cross_product(F,F)(i_3, j_3);
      auto detF =  1./3. * cofF(i_3,j_3) * F(i_3,j_3);
      return {detF, cofF};
    }
 

  inline Tensor2_3d compute_inversion(const Tensor2_3d &F)
  {
    Tensor2_3d iF;
    auto [J, cofF] = compute_determinant_and_cofactor(F);
    iF(i_3,j_3) = 1. / J * cofF(j_3,i_3); 
    return iF;
  }


  inline double compute_I1_C_from_F(const Tensor2_3d &F)
  {
  
    return ( F(k_3,0) * F(k_3,0) + F(k_3,1) * F(k_3,1) + F(k_3,2) * F(k_3,2) );
  }

  inline double compute_I2_C_from_F(const Tensor2_3d &F)
  {
    return ( 0.5 * ( (F(k_3,0) * F(k_3,0) + F(k_3,1) * F(k_3,1) + F(k_3,2) * F(k_3,2)) * (F(k_3,0) * F(k_3,0) + F(k_3,1) * F(k_3,1) + F(k_3,2) * F(k_3,2)) - F(k_3,0) * F(k_3,l_3) * F(m_3,l_3) * F(m_3,0) - F(k_3,1) * F(k_3,l_3) * F(m_3,l_3) * F(m_3,1) - F(k_3,2) * F(k_3,l_3) * F(m_3,l_3) * F(m_3,2) ) );
  }



  inline double compute_I3_C_from_F(const Tensor2_3d &F)
  {
    double J = compute_determinant(F);
    return (J * J);
  }


  inline double compute_I1_C_from_C(const Tensor2_3d &C)
  {   
    return C(i_3,i_3);
  }

  inline double compute_I2_C_from_C(const Tensor2_3d &C)
  {
    double trC = C(i_3,i_3);
    return ( 0.5 * ( trC * trC - ( C(k_3,0) * C(k_3,0) + C(k_3,1) * C(k_3,1) + C(k_3,2) * C(k_3,2) ) ) );
  }

  inline double compute_I3_C_from_C(const Tensor2_3d &C)
  {
    double J2 = compute_determinant(C);
    return J2;
  }
  


  inline Tensor2_3d compute_dI1_C_dF(const Tensor2_3d &F)
  {
    Tensor2_3d dI1_C_dF;
    dI1_C_dF(i_3,j_3) = 2. * F(i_3,j_3); 
    return dI1_C_dF;
  }
 

  inline Tensor2_3d compute_dI2_C_dF(const Tensor2_3d &F)
  {
    Tensor2_3d dI2_C_dF;
    dI2_C_dF(i_3,j_3)=  2. * (F(k_3,0) * F(k_3,0) + F(k_3,1) * F(k_3,1) + F(k_3,2) * F(k_3,2)) * F(i_3,j_3) - 2. * F(i_3,m_3) * F(k_3,m_3) * F(k_3,j_3);
    return dI2_C_dF; 
  }


  

  inline Tensor2_3d compute_dI3_C_dF(const Tensor2_3d &F)
  {
    auto [detF, cofF] = compute_determinant_and_cofactor(F);
    Tensor2_3d dI3_C_dF;
    dI3_C_dF(i_3,j_3) = 2. * detF * cofF(i_3,j_3);
    return dI3_C_dF;
  }


  inline Tensor2_3d compute_dJ_dF(const Tensor2_3d &F)
  {
    return compute_cofactor(F);
  }



  inline Tensor4_3d compute_d2I1_C_dF2(const Tensor2_3d &F)
  {
    Tensor2_3d I(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);
    Tensor4_3d d2I1_C_dF2;
    d2I1_C_dF2(i_3, j_3, k_3, l_3)= 2. *  I( i_3, k_3 ) * I( j_3, l_3 );
    return d2I1_C_dF2;
  
  }


  inline Tensor4_3d compute_d2I2_C_dF2(const Tensor2_3d &F)
  {
    Tensor2_3d I(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);
    Tensor4_3d d2I2_C_dF2;
    d2I2_C_dF2(i_3,j_3,k_3,l_3) = 4. * F( i_3, j_3 ) * F( k_3, l_3 ) - 2. * F( i_3, l_3 ) * F( k_3, j_3 ) +  2. * ( F(m_3,0) * F( m_3, 0 ) + F( m_3, 1 ) * F( m_3, 1 ) + F( m_3, 2 ) * F( m_3, 2 ) ) * I( i_3, k_3 ) * I( j_3, l_3 ) - 2. * F( m_3, l_3 ) * F( m_3, j_3 ) * I( i_3, k_3 ) - 2. * F( i_3, m_3) * F( k_3 , m_3 ) * I( j_3, l_3 );
    return d2I2_C_dF2;
  }




  inline Tensor4_3d compute_d2I3_C_dF2(const Tensor2_3d &F)
  {
    auto [detF, cofF] = compute_determinant_and_cofactor(F);
    Tensor4_3d d2I3_C_dF2;
    d2I3_C_dF2(i_3,j_3,k_3,l_3) = 2. * cofF(i_3,j_3) * cofF(k_3,l_3) + 2 * detF * compute_tensor_cross_product(F)(i_3,j_3,k_3,l_3);
    return d2I3_C_dF2;
  }


  inline double compute_I1_Cdev_from_F(const Tensor2_3d &F)
  {
    return pow(compute_determinant(F), -2./3.) * compute_I1_C_from_F(F);
  }

  inline double compute_I2_Cdev_from_F(const Tensor2_3d &F)
  {
    return pow(compute_determinant(F), -4./3.) * compute_I2_C_from_F(F);
  }

  inline Tensor2_3d compute_dI1_Cdev_dF(const Tensor2_3d &F)
  {
    Tensor2_3d dI1_Cdev_dF;
    auto [J, cofF] = compute_determinant_and_cofactor(F);
    dI1_Cdev_dF(i_3,j_3) = 2. * pow(J, -2./3.) * ( F(i_3,j_3) - compute_I1_C_from_F(F)/ (3. * J) * cofF(i_3,j_3) );
    return dI1_Cdev_dF;
  }

  inline Tensor2_3d compute_dI2_Cdev_dF(const Tensor2_3d &F)
  {
    Tensor2_3d iF, dI2_Cdev_dF, C;
    auto [J, cofF] = compute_determinant_and_cofactor(F);
    C(i_3,j_3) = F(k_3,i_3)*F(k_3,j_3);
    auto I2 = compute_I2_C_from_C(C);

    dI2_Cdev_dF(i_3,j_3) = 2. * pow(compute_determinant(F), -4./3.) * (C(k_3,k_3) * F(i_3,j_3) - F(i_3,k_3) * C(k_3,j_3) - 2. * I2/(3. * J) * cofF(i_3, j_3)  );
    return dI2_Cdev_dF;
  }

 
  inline Tensor4_3d compute_d2I1_Cdev_dF2(const Tensor2_3d &F)
  {
    Tensor2_3d I(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);
    auto I1 = compute_I1_C_from_F(F);
    auto iF = compute_inversion(F);
    Tensor4_3d d2I1_Cdev_dF2;
    d2I1_Cdev_dF2(i_3,j_3,k_3,l_3) = 2./3. * pow(compute_determinant(F), -2./3.) * ( 3. * I( i_3, k_3 ) * I( j_3, l_3 ) + I1 * iF( j_3, k_3 )* iF( l_3, i_3 ) + 2./3. * I1 * iF( j_3, i_3 )* iF( l_3, k_3 )- 2. * iF( l_3, k_3 )* F( i_3, j_3 )- 2. * iF( j_3, i_3 )* F( k_3, l_3 ));
    return d2I1_Cdev_dF2;
  }

  inline Tensor4_3d compute_d2I2_Cdev_dF2(const Tensor2_3d &F) 
  {
    Tensor2_3d I(1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0);
    auto [J, cofF] = compute_determinant_and_cofactor(F);
    auto J_43 = pow(J, -4./3.);
    Tensor2_3d iF, C, B, FC;
    iF(i_3,j_3) = cofF(j_3, i_3)/J;
    C(i_3,j_3) = F(k_3, i_3) * F(k_3, j_3);
    B(i_3,j_3) = F(i_3, k_3) * F(j_3, k_3);
    FC(i_3,j_3) = F(i_3, k_3) * C(k_3, j_3);
    auto I1 = compute_I1_C_from_C(C);
    auto I2 = compute_I2_C_from_C(C);
    Tensor4_3d d2I2_Cdev_dF2;
    d2I2_Cdev_dF2(i_3,j_3,k_3,l_3) = 2. * J_43 * ( I1 * I( i_3, k_3 ) * I( j_3, l_3 ) + 2. * F( i_3, j_3 ) * F( k_3, l_3 ) - 4./3. * I1 * F( i_3, j_3 )  * iF( l_3, k_3 ) + 8./9. * I2 * iF( j_3, i_3 )  * iF( l_3, k_3 ) - 4./3. * I1 * iF( j_3, i_3 )  * F( k_3, l_3 ) + 4. / 3. * FC(k_3,l_3) *  iF(j_3,i_3) + 2./3. * I2 * iF( j_3, k_3 )  * iF( l_3, i_3 ) + 4./3. * FC( i_3, j_3 ) * iF( l_3, k_3 ) - C( l_3, j_3 ) * I(i_3, k_3) - F( i_3, l_3 ) * F( k_3, j_3 ) - B( i_3, k_3 ) * I(l_3, j_3) );
    return  d2I2_Cdev_dF2;
  }



  inline Tensor4_3d compute_d2J_dF2(const Tensor2_3d &F)
  {
    return compute_tensor_cross_product(F);
  }
  

} // end namespace oofem


