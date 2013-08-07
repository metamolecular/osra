/******************************************************************************
 OSRA: Optical Structure Recognition Application

 Created by Igor Filippov, 2007-2013 (igor.v.filippov@gmail.com)

 This program is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
 version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with
 this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
 St, Fifth Floor, Boston, MA 02110-1301, USA
 *****************************************************************************/

// Header: osra_fragments.h
//
// Declares operations on molecular fragments
//
#ifndef OSRA_FRAGMENTS_H
#define OSRA_FRAGMENTS_H

#include "osra.h"

//struct: fragment_s
// used by <osra.cpp::populate_fragments()> to split chemical structure into unconnected molecules.
struct fragment_s
{
  //int: x1,y1,x2,y2
  //top left and bottom right coordinates of the fragment
  int x1, y1, x2, y2;
  //array: atom
  //vector of atom indices for atoms in a molecule of this fragment
  vector<int> atom;
};
//typedef: fragment_t
//defines fragment_t type based on fragment_s struct
typedef struct fragment_s fragment_t;

//
// Section: Functions
//

// Function: find_fragments()
//
// Find disjointed  fragments in a molecule
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
//
// Returns:
// vector of vectors of atom id's which belong to different fragments
vector<vector<int> > find_fragments(const vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom);

// Function: reconnect_fragments()
//
// Reconnecting atoms from different fragments if they are less than 1.1 avg bond length apart
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// avg - average bond length
//
// Returns:
// New number of bonds
int reconnect_fragments(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double avg);

// Function: populate_fragments()
//
// Transforms vector of vectors of atom ids into a vector of fragments
//
// Parameters:
// frags - vector of vector of atom ids
// atom - vector of atoms
//
// Returns:
// vector of fragments
vector<fragment_t> populate_fragments(const vector<vector<int> > &frags, const vector<atom_t> &atom);

// Function: comp_fragments()
//
// Comparison function used for sorting fragments according to their positions in the picture: top-down, left to right
//
// Parameters:
// aa, bb - two fragments to compare
//
// Returns:
// True if fragment aa is higher or to the left of fragment bb.
// False otherwise.
bool comp_fragments(const fragment_t &aa, const fragment_t &bb);

#endif // OSRA_FRAGMENTS_H
