// $Id: ESMCI_MeshSkin.h,v 1.5 2011/01/05 20:05:44 svasquez Exp $
// Earth System Modeling Framework
// Copyright 2002-2011, University Corporation for Atmospheric Research, 
// Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
// Laboratory, University of Michigan, National Centers for Environmental 
// Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
// NASA Goddard Space Flight Center.
// Licensed under the University of Illinois-NCSA License.

//
//-----------------------------------------------------------------------------
#ifndef ESMCI_MeshSkin_h
#define ESMCI_MeshSkin_h

#include "ESMCI_Mesh.h"

namespace ESMCI {

// Add sidesets for exterior blocks and for interior block boundaries
// (Locally) mark nodes as boundary nodes.  If add_fadces=true, create the faces.
void Skin(Mesh &mesh);

// For a parallel mesh, resolve sides on a parallel boundary, i.e. delete sides
// that have been created on a parallel boundary because the local proc thinks this
// is an exposed boundary.  Also, deimprint the exposedboundary from the corresponding
// nodes.
void ResolveParSkin(Mesh &mesh);

} //namespace

#endif
