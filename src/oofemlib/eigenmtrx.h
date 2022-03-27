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



#ifndef eigenmtrx_h
#define eigenmtrx_h

#include "sparsemtrx.h"
#include "intarray.h"

#include <Eigen/Sparse>

#define _IFT_EigenMtrx_Name "EigenMtrx"

namespace oofem {
/**
 * Implementation of sparse matrix stored using Eigen library.
 */
class OOFEM_EXPORT EigenMtrx : public SparseMtrx
{
protected:
    Eigen::SparseMatrix<double> EigMat;
    //std::vector<Eigen::Triplet<double> > triplets; // Allocate vector of triplets

    //FloatArray val; // data values (nz_ elements)
    //IntArray rowind; // row_ind (nz_ elements)
    //IntArray colptr; // col_ptr (dim_[1]+1 elements)

    //int base; // index base: offset of first element
    //int nz; // number of nonzeros

public:
    /** Constructor. Before any operation an internal profile must be built.
     * @see buildInternalStructure
     */
    EigenMtrx( int n=0 );

    /// Destructor
    virtual ~EigenMtrx() {}

    // Overloaded methods:
    /*std::unique_ptr<SparseMtrx> clone() const override;
    void times( const FloatArray &x, FloatArray &answer ) const override;
    void timesT( const FloatArray &x, FloatArray &answer ) const override;
    void times( double x ) override;*/
    int buildInternalStructure( EngngModel *, int, const UnknownNumberingScheme &s ) override;
    int assemble( const IntArray &loc, const FloatMatrix &mat ) override;
    int assemble( const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat ) override;
    bool canBeFactorized() const override { return false; }
    void zero() override;
    double &at( int i, int j ) override;
    double at( int i, int j ) const override;
    /*void toFloatMatrix( FloatMatrix &answer ) const override;
    void printYourself() const override;*/
    const char *giveClassName() const override { return "EigenMtrx"; }
    SparseMtrxType giveType() const override { return SMT_EigenMtrx; }
    bool isAsymmetric() const override { return true; }

    Eigen::SparseMatrix<double> giveMatrix();



};
} // end namespace oofem
#endif // eigenmtrx_h

