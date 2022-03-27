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

#include "eigensolver.h"
#include "eigenmtrx.h"
//#include <iostream>

#include "compcol.h"
#include "symcompcol.h"
#include "engngm.h"
#include "floatarray.h"
#include "verbose.h"
#include "timer.h"
#include "error.h"
#include "classfactory.h"
//
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

namespace oofem {
REGISTER_SparseLinSolver( EigenSolver, ST_Eigen );

EigenSolver ::EigenSolver( Domain *d, EngngModel *m ) :
    SparseLinearSystemNM( d, m ) 
{
}

EigenSolver ::~EigenSolver() {}

IRResultType
EigenSolver :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // required by IR_GIVE_FIELD macro
    this->method = "llt";
    IR_GIVE_OPTIONAL_FIELD( ir, method, "eigenmethod" );
    return IRRT_OK;

}
  

NM_Status EigenSolver ::solve( SparseMtrx &A, FloatArray &b, FloatArray &x )
{
    int neqs = b.giveSize(); // Number of equations

    EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &A );

    Eigen::SparseMatrix<double> A_eig = Ae->giveMatrix();

    // Construct right hand side vetor
    Eigen::VectorXd b_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( b.givePointer(), neqs );

    Timer timer;
    timer.startTimer();


    Eigen::VectorXd x_eig; // Allocate vector of RHS

    // Create factorization
    if ( this->method.compare("LLT" ) == 0 ) {

        Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > A_factorization( A_eig );       
        x_eig = A_factorization.solve( b_eig ); // Solve the system

    } else if ( this->method.compare("LDLT" ) == 0 ) {

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > A_factorization( A_eig );     
        x_eig = A_factorization.solve( b_eig ); // Solve the system

    } else if ( this->method.compare("LU" ) == 0 ) {

        Eigen::SparseLU<Eigen::SparseMatrix<double> > A_factorization( A_eig );      
        x_eig = A_factorization.solve( b_eig ); // Solve the system

    } else if ( this->method.compare("QR" ) == 0 ) {

        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > A_factorization( A_eig );
        x_eig = A_factorization.solve( b_eig ); // Solve the system
    }else if ( this->method.compare("LDLT_Dcontrol" ) == 0 ) {

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > A_factorization( A_eig );
	A_factorization.updateD();
        x_eig = A_factorization.solve( b_eig ); // Solve the system

    } else {
      OOFEM_ERROR("Unknown solver type");
    }

    // Copy/move values to FloatArray x
    x = FloatArray( x_eig.begin(), x_eig.end() );

    timer.stopTimer();
    OOFEM_LOG_INFO( "EigenSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );

    NM_Status s = NM_Success;
    return s;

}

} // end namespace oofem

