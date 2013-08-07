/******************************************************************************
 OSRA: Optical Structure Recognition

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

#ifndef OSRA_OPENBABEL_H
#define OSRA_OPENBABEL_H

#include <string> // std::string
#include <map> // std::map
#include <vector> // std::vector

#include "osra.h"
#include "osra_segment.h"

using namespace std;

// Header: osra_openbabel.h
//
// Defines types and functions for OSRA OpenBabel module.
//


//struct: molecule_statistics_s
//      contains the statistical information about molecule used for analysis of recognition accuracy
struct molecule_statistics_s
{
  // int: rotors
  //    number of rotors in molecule
  int rotors;
  // int: num_fragments
  //    number of contiguous fragments in molecule
  int fragments;
  // int: rings56
  //    accumulated number of 5- and 6- rings in molecule
  int rings56;
  // int: num_atoms
  // number of atoms in molecule
  int num_atoms;
// int: num_bonds
  // number of bonds in molecule
  int num_bonds;
};

// typedef: molecule_statistics_t
//      defines molecule_statistics_t type based on molecule_statistics_s struct
typedef struct molecule_statistics_s molecule_statistics_t;

//
// Section: Functions
//

// Function: osra_openbabel_init()
//
// Performs OpenBabel library engine sanity check. Should be called at e.g. program startup.
//
// Returns:
//      non-zero value in case of error
int osra_openbabel_init();

// Function: calculate_molecule_statistics()
//
// Converts vectors of atoms and bonds into a molecular object and calculates the molecule statistics.
// Note: this function changes the atoms!
//
// Parameters:
//      atom - vector of <atom_s> atoms
//      bond - vector of <bond_s> bonds
//      n_bond - total number of bonds
//      avg_bond_length - average bond length as measured from the image (to be included into output if provided)
//      superatom - dictionary of superatom labels mapped to SMILES
//      verbose - print debug info
//
// Returns:
//      calculated molecule statistics
molecule_statistics_t calculate_molecule_statistics(vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond,
    double avg_bond_length, const map<string, string> &superatom, bool verbose);

// Function: get_formatted_structure()
//
// Converts vectors of atoms and bonds into a molecular object and encodes the molecular into a text presentation (SMILES, MOL file, ...),
// specified by given format.
//
// Parameters:
//      atom - vector of <atom_s> atoms
//      bond - vector of <bond_s> bonds
//      n_bond - total number of bonds
//      format - output format for molecular representation - i.e. SMI, SDF
//      embedded_format - output format to be embedded into SDF (is only valid if output format is SDF); the only embedded formats supported now are "inchi", "smi", and "can"
//      molecule_statistics - the molecule statistics (returned to the caller)
//      confidence - confidence score (returned to the caller)
//      show_confidence - toggles confidence score inclusion into output
//      avg_bond_length - average bond length as measured from the image
//      scaled_avg_bond_length - average bond length scaled to the original resolution of the image
//      show_avg_bond_length - toggles average bond length inclusion into output
//      resolution - resolution at which image is being processed in DPI (to be included into output if provided)
//      page - page number (to be included into output if provided)
//      surrounding_box - the coordinates of surrounding image box that contains the structure (to be included into output if provided)
//      superatom - dictionary of superatom labels mapped to SMILES
//      verbose - print debug info
//
//  Returns:
//      string containing SMILES, SDF or other representation of the molecule
const string get_formatted_structure(vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, const string &format, const string &second_format, molecule_statistics_t &molecule_statistics,
                                     double &confidence, bool show_confidence, double avg_bond_length, double scaled_avg_bond_length, bool show_avg_bond_length, const int * const resolution,
                                     const int * const page, const box_t * const surrounding_box, const map<string, string> &superatom, int n_letters, bool show_learning, int resolution_iteration, bool verbose);

#endif
