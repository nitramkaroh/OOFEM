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


#include <ostream>

namespace oofem {

/**
* Implementation of the second order tensor class
* @author Martin Horak, add yourself
*/
  
template <typename Type = double, int dimension1 = 3, int dimension2 = 3 >
class OOFEM_EXPORT Tensor2
{
 protected:
 /// Values of matrix stored column wise.
 std::array< Type, Dimension_1*Dimension_2 > values;
 int rank = 2;
 
 public:
 /// @name Iterator for for-each loops:
 //@{
 auto begin() { return this->values.begin(); }
 auto end() { return this->values.end(); }
 auto begin() const { return this->values.begin(); }
 auto end() const { return this->values.end(); }
 //@}
 constexpr inline Tensor2() = default;
/**
     * Constructor (values are specified column-wise)
     * @note The syntax {{x,y,z},{...}} can be achieved by nested initializer_list, but 
     */
    template< typename... V, class = typename std::enable_if_t<sizeof...(V) == Dimension_1*Dimension_2> >
    Tensor2(V... x) noexcept : values{x...} { }


    /**
     * Empty ctor, initializes to zero.
     */
     Tensor2() noexcept : values{} { }
    /// Copy constructor.
    Tensor2(const Tensor2<type, Dimension_1,Dimension_2> &ten2) noexcept : values(ten2.values) {}
    /// FloatMatrix conversion constructor.
    Tensor2(const FloatMatrix &mat)
    {
#ifndef NDEBUG
        if ( mat.giveNumberOfRows() != Dimension_1 || mat.giveNumberOfColumns() != Dimension_2 ) {
            throw std::out_of_range("Can't convert dynamic float matrix of size " + 
                std::to_string(mat.giveNumberOfRows()) + "x" + std::to_string(mat.giveNumberOfRows()) + 
                " to fixed size " + std::to_string(N) + "x" + std::to_string(M));
        }
#endif
        std::copy_n(mat.begin(), Dimension_1*Dimension_2, values.begin());
    }
    Tensor2(Tensor2<Dimension_1> const (&x)[Dimension_2]) noexcept
    {
        for (std::size_t i = 0; i < Dimension_1; ++i) {
            for (std::size_t j = 0; j < Dimension_2; ++j) {
                (*this)(i,j) = x[j][i];
            }
        }
    }
 
    /// Assignment operator
    Tensor2 &operator=(const Tensor2<Dimension_1,Dimension_2> &ten)
    {
        values = ten.values;
        return * this;
    }

 /**
  * Symmetrizes and stores a matrix in Voigt form:
  * x_11, x_22, x_33, x_23, x_13, x_12, x_32, x_31, x_21
  */
 inline Tensor2<3,3> from_voigt_form(const FloatArrayF<9> &v)
 {
   return {
     v[0], v[8], v[7], 
     v[5], v[1], v[6], 
     v[4], v[3], v[2]
   }; 
 }


 double &operator()(int i, int j)
 {
   return values [ j * Dimension_1 + i ];
 }
 
};
 
}


