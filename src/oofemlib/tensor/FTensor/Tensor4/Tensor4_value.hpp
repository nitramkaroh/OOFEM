/* A general version, not for pointers. */

#pragma once

#include "Tensor4_contracted.hpp"

#include <ostream>
#include "floatmatrixf.h"
#ifdef FTENSOR_DEBUG
#include <sstream>
#include <stdexcept>
#endif

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3>
  class Tensor4
  {
  protected:
    T data[Tensor_Dim0][Tensor_Dim1][Tensor_Dim2][Tensor_Dim3];

  public:
    /* Initialization operators */
    template <class... U> Tensor4(U... d) : data{d...}
    {
      static_assert(sizeof...(d) == sizeof(data) / sizeof(T),
                    "Incorrect number of Arguments. Constructor should "
                    "initialize the entire Tensor");
    }

    Tensor4() {}
    /* There are two operator(int,int,int,int)'s, one for non-consts
       that lets you change the value, and one for consts that
       doesn't. */
    T &operator()(const int N1, const int N2, const int N3, const int N4)
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
         || N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor4<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << "," << Tensor_Dim2 << "," << Tensor_Dim3 << ">.operator(" << N1
            << "," << N2 << "," << N3 << "," << N4 << ") const" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N1][N2][N3][N4];
    }

    T operator()(const int N1, const int N2, const int N3, const int N4) const
    {
#ifdef FTENSOR_DEBUG
      if(N1 >= Tensor_Dim0 || N1 < 0 || N2 >= Tensor_Dim1 || N2 < 0
         || N3 >= Tensor_Dim2 || N3 < 0 || N4 >= Tensor_Dim3 || N4 < 0)
        {
          std::stringstream s;
          s << "Bad index in Tensor4<T," << Tensor_Dim0 << "," << Tensor_Dim1
            << "," << Tensor_Dim2 << "," << Tensor_Dim3 << ">.operator(" << N1
            << "," << N2 << "," << N3 << "," << N4 << ")" << std::endl;
          throw std::out_of_range(s.str());
        }
#endif
      return data[N1][N2][N3][N4];
    }

    /* These operator()'s are the first part in constructing template
       expressions.  They can be used to slice off lower dimensional
       parts. */

    template <char i, char j, char k, char l, int Dim0, int Dim1, int Dim2,
              int Dim3>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2
       && Tensor_Dim3 >= Dim3),
      Tensor4_Expr<Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>,
                   T, Dim0, Dim1, Dim2, Dim3, i, j, k, l>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
               const Index<k, Dim2>, const Index<l, Dim3>)
    {
      return Tensor4_Expr<
        Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>, T,
        Dim0, Dim1, Dim2, Dim3, i, j, k, l>(*this);
    };

    template <char i, char j, char k, char l, int Dim0, int Dim1, int Dim2,
              int Dim3>
    typename std::enable_if<
      (Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1 && Tensor_Dim2 >= Dim2
       && Tensor_Dim3 >= Dim3),
      Tensor4_Expr<
        const Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>,
        T, Dim0, Dim1, Dim2, Dim3, i, j, k, l>>::type
    operator()(const Index<i, Dim0>, const Index<j, Dim1>,
               const Index<k, Dim2>, const Index<l, Dim3>) const
    {
      return Tensor4_Expr<
        const Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>,
        T, Dim0, Dim1, Dim2, Dim3, i, j, k, l>(*this);
    };

    /* These operators are for internal contractions, resulting in a Tensor2.
     * For example something like A(k,i,k,j) */

    // (i,j,k,k)
    template <char i, char j, char k, int Dim0, int Dim1, int Dim23>
    auto operator()(const Index<i, Dim0>, const Index<j, Dim1>,
                    const Index<k, Dim23>, const Index<k, Dim23>) const
    {
      static_assert(Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim1
                      && Tensor_Dim2 >= Dim23 && Tensor_Dim3 >= Dim23,
                    "Incompatible indices");

      using TensorExpr = Tensor4_contracted_23<
        Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>, T,
        Dim23>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim1, i, j>(TensorExpr(*this));
    };

    // (i,j,k,j)
    template <char i, char j, char k, int Dim0, int Dim13, int Dim2>
    auto operator()(const Index<i, Dim0>, const Index<j, Dim13>,
                    const Index<k, Dim2>, const Index<j, Dim13>) const
    {
      static_assert(Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim13
                      && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim13,
                    "Incompatible indices");

      using TensorExpr = Tensor4_contracted_13<
        Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>, T,
        Dim13>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim2, i, k>(TensorExpr(*this));
    };

    // (i,j,j,l)
    template <char i, char j, char l, int Dim0, int Dim12, int Dim3>
    auto operator()(const Index<i, Dim0>, const Index<j, Dim12>,
                    const Index<j, Dim12>, const Index<l, Dim3>) const
    {
      static_assert(Tensor_Dim0 >= Dim0 && Tensor_Dim1 >= Dim12
                      && Tensor_Dim2 >= Dim12 && Tensor_Dim3 >= Dim3,
                    "Incompatible indices");

      using TensorExpr = Tensor4_contracted_12<
        Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>, T,
        Dim12>;
      return Tensor2_Expr<TensorExpr, T, Dim0, Dim3, i, l>(TensorExpr(*this));
    };

    // (i,j,k,i)
    template <char i, char j, char k, int Dim03, int Dim1, int Dim2>
    auto operator()(const Index<i, Dim03>, const Index<j, Dim1>,
                    const Index<k, Dim2>, const Index<i, Dim03>) const
    {
      static_assert(Tensor_Dim0 >= Dim03 && Tensor_Dim1 >= Dim1
                      && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim03,
                    "Incompatible indices");

      using TensorExpr = Tensor4_contracted_03<
        Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>, T,
        Dim03>;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim2, j, k>(TensorExpr(*this));
    };

    // (i,j,i,l)
    template <char i, char j, char l, int Dim02, int Dim1, int Dim3>
    auto operator()(const Index<i, Dim02>, const Index<j, Dim1>,
                    const Index<i, Dim02>, const Index<l, Dim3>) const
    {
      static_assert(Tensor_Dim0 >= Dim02 && Tensor_Dim1 >= Dim1
                      && Tensor_Dim2 >= Dim02 && Tensor_Dim3 >= Dim3,
                    "Incompatible indices");

      using TensorExpr = Tensor4_contracted_02<
        Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>, T,
        Dim02>;
      return Tensor2_Expr<TensorExpr, T, Dim1, Dim3, j, l>(TensorExpr(*this));
    };

    // (i,i,k,l)
    template <char i, char k, char l, int Dim01, int Dim2, int Dim3>
    auto operator()(const Index<i, Dim01>, const Index<i, Dim01>,
                    const Index<k, Dim2>, const Index<l, Dim3>) const
    {
      static_assert(Tensor_Dim0 >= Dim01 && Tensor_Dim1 >= Dim01
                      && Tensor_Dim2 >= Dim2 && Tensor_Dim3 >= Dim3,
                    "Incompatible indices");

      using TensorExpr = Tensor4_contracted_01<
        Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3>, T,
        Dim01>;
      return Tensor2_Expr<TensorExpr, T, Dim2, Dim3, k, l>(TensorExpr(*this));
    };

    /* This is for expressions where a number is used for one slot, and
       an index for another, yielding a Tensor3_Expr.  The non-const
       versions don't actually create a Tensor3_number_rhs_[0123] object.
       They create a Tensor3_Expr directly, which provides the
       appropriate indexing operators.  The const versions do create a
       Tensor3_number_[0123]. */
    // TODO
  };


  class Tensor4_3d : public Tensor4<double, 3, 3, 3, 3>
  {
  public:
    using Tensor4 <double, 3, 3, 3, 3> :: Tensor4;
    Tensor4_3d(const oofem::FloatMatrixF<9,9> &mat){
      // the first column
      this->data[0][0][0][0] = mat.at(1,1);
      this->data[0][1][0][0] = mat.at(6,1);
      this->data[0][2][0][0] = mat.at(5,1);
      this->data[1][0][0][0] = mat.at(9,1);
      this->data[1][1][0][0] = mat.at(2,1);
      this->data[1][2][0][0] = mat.at(4,1);
      this->data[2][0][0][0] = mat.at(8,1);
      this->data[2][1][0][0] = mat.at(7,1);
      this->data[2][2][0][0] = mat.at(3,1);
      // the second column
      this->data[0][0][1][1] = mat.at(1,2);
      this->data[0][1][1][1] = mat.at(6,2);
      this->data[0][2][1][1] = mat.at(5,2);
      this->data[1][0][1][1] = mat.at(9,2);
      this->data[1][1][1][1] = mat.at(2,2);
      this->data[1][2][1][1] = mat.at(4,2);
      this->data[2][0][1][1] = mat.at(8,2);
      this->data[2][1][1][1] = mat.at(7,2);
      this->data[2][2][1][1] = mat.at(3,2);
      // the third column
      this->data[0][0][2][2] = mat.at(1,3);
      this->data[0][1][2][2] = mat.at(6,3);
      this->data[0][2][2][2] = mat.at(5,3);
      this->data[1][0][2][2] = mat.at(9,3);
      this->data[1][1][2][2] = mat.at(2,3);
      this->data[1][2][2][2] = mat.at(4,3);
      this->data[2][0][2][2] = mat.at(8,3);
      this->data[2][1][2][2] = mat.at(7,3);
      this->data[2][2][2][2] = mat.at(3,3);
      // the fourth column
      this->data[0][0][1][2] = mat.at(1,4);
      this->data[0][1][1][2] = mat.at(6,4);
      this->data[0][2][1][2] = mat.at(5,4);
      this->data[1][0][1][2] = mat.at(9,4);
      this->data[1][1][1][2] = mat.at(2,4);
      this->data[1][2][1][2] = mat.at(4,4);
      this->data[2][0][1][2] = mat.at(8,4);
      this->data[2][1][1][2] = mat.at(7,4);
      this->data[2][2][1][2] = mat.at(3,4);
      //the fifth column
      this->data[0][0][0][2] = mat.at(1,5);
      this->data[0][1][0][2] = mat.at(6,5);
      this->data[0][2][0][2] = mat.at(5,5);
      this->data[1][0][0][2] = mat.at(9,5);
      this->data[1][1][0][2] = mat.at(2,5);
      this->data[1][2][0][2] = mat.at(4,5);
      this->data[2][0][0][2] = mat.at(8,5);
      this->data[2][1][0][2] = mat.at(7,5);
      this->data[2][2][0][2] = mat.at(3,5);
      //the sixth column
      this->data[0][0][0][1] = mat.at(1,6);
      this->data[0][1][0][1] = mat.at(6,6);
      this->data[0][2][0][1] = mat.at(5,6);
      this->data[1][0][0][1] = mat.at(9,6);
      this->data[1][1][0][1] = mat.at(2,6);
      this->data[1][2][0][1] = mat.at(4,6);
      this->data[2][0][0][1] = mat.at(8,6);
      this->data[2][1][0][1] = mat.at(7,6);
      this->data[2][2][0][1] = mat.at(3,6);
      //the seventh column
      this->data[0][0][2][1] = mat.at(1,7);
      this->data[0][1][2][1] = mat.at(6,7);
      this->data[0][2][2][1] = mat.at(5,7);
      this->data[1][0][2][1] = mat.at(9,7);
      this->data[1][1][2][1] = mat.at(2,7);
      this->data[1][2][2][1] = mat.at(4,7);
      this->data[2][0][2][1] = mat.at(8,7);
      this->data[2][1][2][1] = mat.at(7,7);
      this->data[2][2][2][1] = mat.at(3,7);
      //the eight column
      this->data[0][0][2][0] = mat.at(1,8);
      this->data[0][1][2][0] = mat.at(6,8);
      this->data[0][2][2][0] = mat.at(5,8);
      this->data[1][0][2][0] = mat.at(9,8);
      this->data[1][1][2][0] = mat.at(2,8);
      this->data[1][2][2][0] = mat.at(4,8);
      this->data[2][0][2][0] = mat.at(8,8);
      this->data[2][1][2][0] = mat.at(7,8);
      this->data[2][2][2][0] = mat.at(3,8);
      //the nineght column
      this->data[0][0][1][0] = mat.at(1,9);
      this->data[0][1][1][0] = mat.at(6,9);
      this->data[0][2][1][0] = mat.at(5,9);
      this->data[1][0][1][0] = mat.at(9,9);
      this->data[1][1][1][0] = mat.at(2,9);
      this->data[1][2][1][0] = mat.at(4,9);
      this->data[2][0][1][0] = mat.at(8,9);
      this->data[2][1][1][0] = mat.at(7,9);
      this->data[2][2][1][0] = mat.at(3,9);
    }
    

    void __attribute__((noinline)) printYourself() const 
    // Prints the receiver on screen.
    {
      printf("4th order tensor");
      for ( int i = 1; i <= 3; i++ ) {
	for ( int j = 1; j <= 3; j++ ) {
	  for ( int k = 1; k <= 3; k++ ) {
	    for ( int l = 1; l <= 3; l++ ) {
	      printf( "Component %d %d %d %d is: %10.3e \n", i,j,k,l,this->data[i][j][k][l] );
	    }
	  }
	}
      }  
    }
    
    void __attribute__((noinline)) printComponent(int i, int j, int k, int l) const 
    // Prints the receiver on screen.
    {
      printf( "Component %d %d %d %d is: %10.3e\n  ", i,j,k,l,this->data[i][j][k][l] );
    }
    

  };

}

/// JSON compatible output

namespace FTensor
{
  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3>
  std::ostream &Tensor4_0001(
    std::ostream &os,
    const Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3> &t,
    const int &iterator0, const int &iterator1, const int &iterator2)
  {
    os << '[';
    for(int i = 0; i < Tensor_Dim3 - 1; ++i)
      {
        os << t(iterator0, iterator1, iterator2, i);
        os << ',';
      }
    os << t(iterator0, iterator1, iterator2, Tensor_Dim3 - 1);
    os << ']';

    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3>
  std::ostream &Tensor4_0010(
    std::ostream &os,
    const Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3> &t,
    const int &iterator0, const int &iterator1)
  {
    os << '[';
    for(int i = 0; i < Tensor_Dim2 - 1; ++i)
      {
        FTensor::Tensor4_0001(os, t, iterator0, iterator1, i);
        os << ',';
      }
    FTensor::Tensor4_0001(os, t, iterator0, iterator1, Tensor_Dim2 - 1);
    os << ']';

    return os;
  }

  template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
            int Tensor_Dim3>
  std::ostream &Tensor4_0100(
    std::ostream &os,
    const Tensor4<T, Tensor_Dim0, Tensor_Dim1, Tensor_Dim2, Tensor_Dim3> &t,
    const int &iterator0)
  {
    os << '[';
    for(int i = 0; i < Tensor_Dim1 - 1; ++i)
      {
        FTensor::Tensor4_0010(os, t, iterator0, i);
        os << ',';
      }
    FTensor::Tensor4_0010(os, t, iterator0, Tensor_Dim1 - 1);
    os << ']';

    return os;
  }
}

template <class T, int Tensor_Dim0, int Tensor_Dim1, int Tensor_Dim2,
          int Tensor_Dim3>
std::ostream &operator<<(std::ostream &os,
                         const FTensor::Tensor4<T, Tensor_Dim0, Tensor_Dim1,
                                                Tensor_Dim2, Tensor_Dim3> &t)
{
  os << '[';
  for(int i = 0; i < Tensor_Dim0 - 1; ++i)
    {
      FTensor::Tensor4_0100(os, t, i);
      os << ',';
    }
  FTensor::Tensor4_0100(os, t, Tensor_Dim0 - 1);
  os << ']';

  return os;
}
