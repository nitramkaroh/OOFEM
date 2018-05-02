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

#ifndef podvtkxmlexportmodule_h
#define podvtkxmlexportmodule_h

#include "vtkxmlexportmodule.h"


///@name Input fields for VTK XML export module
//@{
#define _IFT_PODVTKXMLExportModule_Name "podvtkxml"
#define _IFT_PODVTKXMLExportModule_primvars "primvars"
#define _IFT_PODVTKXMLExportModule_exportFileName "exportfilename"

//@}

namespace oofem {
class Node;

/**
 * Represents VTK (Visualization Toolkit) export module. It uses VTK (.vtu) file format, Unstructured grid dataset.
 * The export of data is done on Region By Region basis, possibly taking care about possible nonsmooth character of
 * some internal variables at region boundaries.
 * Each region is usually exported as a single piece. When region contains composite cells, these are assumed to be
 * exported in individual subsequent pieces after the default one for the particular region.
 */
class OOFEM_EXPORT PODVTKXMLExportModule : public VTKXMLExportModule
{
protected:
  std :: string exportFileName;

public:
    /// Constructor. Creates empty Output Manager. By default all components are selected.
    PODVTKXMLExportModule(int n, EngngModel * e);
    /// Destructor
    virtual ~PODVTKXMLExportModule();

    virtual void doOutput(TimeStep *tStep, int nSnapshot);
    virtual const char *giveClassName() const { return "PODVTKXMLExportModule"; }
    virtual void  writeVTKPiece(VTKPiece &vtkPiece, TimeStep *tStep);
    IRResultType initializeFrom(InputRecord *ir);
    virtual void initialize();
    




protected:

    /// Returns the filename for the given time step.
    std :: string giveOutputFileName(int iSnapshot);

    /// Returns the output stream for given solution step.
    FILE *giveOutputStream(int iSnapshot);

    void doOutputSnapshot(TimeStep *tStep, int iSnapshot);
    void exportPrimaryVars(VTKPiece &piece, IntArray &mapG2L, IntArray &mapL2G, int region, TimeStep *tStep, int iSnapshot);
    void getNodalVariableFromPrimaryField(FloatArray &answer, DofManager *dman, TimeStep *tStep, UnknownType type, int ireg, int iSnapshot);
    void setupVTKPiece(VTKPiece &vtkPiece,TimeStep *tStep, int region, int iSnapshot);
};





} // end namespace oofem
#endif // podvtkxmlexportmodule_h
