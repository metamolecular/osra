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

// Header: osra_structure.h
//
// Declares main structure recognition (molecular atoms and bonds)
// recognition routines
//

#include "osra.h"
#include "osra_labels.h"



//struct: dash_s
// used to identify dashed bonds in <osra.cpp::find_dashed_bonds()> and small bonds in <osra.cpp::find_small_bonds()>
struct dash_s
{
  //double: x,y
  //coordinates
  double x, y;
  //bool: free
  // is this dash available for a perspective dashed bond?
  bool free;
  //pointer: curve
  // pointer to the curve found by Potrace
  const potrace_path_t *curve;
  //int: area
  // area occupied by the dash
  int area;
};
//typedef: dash_t
//defines dash_t type based on dash_s struct
typedef struct dash_s dash_t;

//
// Section: Functions
//

// Function: remove_disconnected_atoms()
//
// Removes atoms not connected to any bond
//
// Parameters:
// atom - vector of atoms
// bond - vector of bonds
// n_atom - number of atoms
// n_bond - number of bonds
void remove_disconnected_atoms(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond);

// Function: remove_zero_bonds()
//
// Removes bonds which either coincide with an existing bond, or have beginning and end atoms the same, or
// have one of the atoms non-existing
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
void remove_zero_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom);

// Function: collapse_doubleup_bonds()
//
// if a bond coincide with another, remove it and increment the bond type of the existing bond
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
void collapse_doubleup_bonds(vector<bond_t> &bond, int n_bond);

// Function: skeletize()
//
// Collapses two sides of the same bond into a single bond object
//
// Parameters:
// atom - vector of atoms
// bond - vector of bonds
// n_bond - number of bonds
// image - original image
// threshold - black-white binarization threshold
// bgColor - background color
// dist - distance below which it's considered the same bond no matter whether the  other conditions are met or not
// avg - average bond length
//
// Returns:
// Average bond thickness
double skeletize(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, const Image &image, double threshold,const ColorGray &bgColor, double dist, double avg);

// Function: dist_double_bonds()
//
// Estimates maximum distance between double bonds
//
// Parameters:
// atom - vector of atoms
// bond - vector of bonds
// n_bond - number of bonds
// avg - average bond length
//
// Returns:
// Estimated distance between double bonds
double dist_double_bonds(const vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, double avg);

// Function: double_triple_bonds()
//
// Detects double and triple bonds from single collinear bonds
//
// Parameters:
// atom - vector of atoms
// bond - vector of bonds
// n_bond - number of bonds
// avg - average bond length
// n_atom - number of atoms
// max_dist_double_bond - maximum distance between double bonds
//
// Returns:
// New value for n_bond
int double_triple_bonds(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, double avg, int &n_atom,double max_dist_double_bond);

// Function: extend_terminal_bond_to_label()
//
// Attempt to connect a terminal bond to an atomic label
//
// Parameters:
// atom - vector of atoms
// letters - vector of characters
// n_letters - number of characters
// bond - vector of bonds
// n_bond - number of bonds
// label - vector of atomic labels
// n_label - number of atomic labels
// avg - average bond length
// maxh - maximum bond thickness
// max_dist_double_bond - maximum distance between double bonds
void extend_terminal_bond_to_label(vector<atom_t> &atom, const vector<letters_t> &letters, int n_letters, const vector<bond_t> &bond, int n_bond, const vector<label_t> &label, int n_label, double avg, double maxh,
                                   double max_dist_double_bond);

// Function: extend_terminal_bond_to_bonds()
//
// Attempts to reconnect a terminal bond to other bonds
//
// Parameters:
// atom - vector of atoms
// bond - vector of bonds
// n_bond - number of bonds
// avg - average bond length
// maxh - maximum bond thickness
// max_dist_double_bond - maximum distance between double bonds
void extend_terminal_bond_to_bonds(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, double avg, double maxh,double max_dist_double_bond);

// Function: assign_charge()
//
// Assignes charges to atoms where recognized. <fix_atom_name()> is also called at the end.
//
// Parameters:
// atom - vector of atoms
// bond - vector of bonds
// n_atom - number of atoms
// n_bond - number of bonds
// fix - spelling correction map loaded from an external file
// superatom - superatom dictionary map loaded from an external file
// debug - a flag for debug output
void assign_charge(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond, const map<string, string> &fix,const map<string, string> &superatom, bool debug);

// Function: find_atoms()
//
// Detects atoms and bonds in a vectorized image
//
// Parameters:
// p - output of Potrace vectorization routines
// atom - vector of atoms
// bond - vector of bonds
// n_bond - pointer to the number of bonds
// width - image width
// height - image height
//
// Returns:
// number of atoms
int find_atoms(const potrace_path_t *p, vector<atom_t> &atom, vector<bond_t> &bond, int *n_bond,int width, int height);

// Function: find_dashed_bonds()
//
// Detects dashed (hash) stereo bonds
//
// Parameters:
// p - output of Potrace vectorization routines
// atom - vector of atoms
// bond - vector of bonds
// n_atom - number of atoms
// n_bond - pointer to the number of bonds
// max - maximum area for a dash
// avg - average bond length
// img - original image
// bg - background color
// THRESHOLD - black-white binarization threshold
// thick - flag set if the image was subject to thinning
// dist - maximum dashed bond length
int find_dashed_bonds(const potrace_path_t *p, vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int *n_bond,int max, double avg, const Image &img,
                      const ColorGray &bg, double THRESHOLD, bool thick, double dist, vector<letters_t> &letters);

// Function: find_small_bonds()
//
// Bonds which are small in size are flagged as such
//
// Parameters:
// p - Potrace vectorization output
// atom - vector of atoms
// bond - vector of bonds
// n_atom - number of atoms
// n_bond - pointer to the number of bonds
// max_area - maximum area of the bond segment
// Small - maximum area of the bond segment which allows not to take into account thickness
// thickness - maximum thickness for a small bond
//
// Returns:
// new value of n_atom
int find_small_bonds(const potrace_path_t *p, vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int *n_bond,
                     double max_area, double Small, double thickness);

// Function: resolve_bridge_bonds()
//
// Detects bridge bonds - bonds which are visually intersecting but have no atom at the point of intersection
//
// Parameters:
// atom - vector of atoms
// n_atom - number of atoms
// bond - vector of bonds
// n_bond - number of bonds
// thickness - bond thickness
// avg_bond_length - average bond length
// supertatom - superatom dictionary map
// verbose - print debug info
//
// Returns:
// The number of fragments
int resolve_bridge_bonds(vector<atom_t> &atom, int n_atom, vector<bond_t> &bond, int n_bond, double thickness,
                         double avg_bond_length, const map<string, string> &superatom, bool verbose);

// Function: collapse_atoms()
//
// Atoms found within a threshold distance from each other are collapsed into one
//
// Parameters:
// atom - vector of atoms
// bond - vector of bonds
// n_atom - number of atoms
// n_bond - number of bonds
// dist - threshold distance between atoms
void collapse_atoms(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond, double dist);

// Function: collapse_bonds()
//
// Bond ends with bond length smaller than a threshold are collapsed into one atom
//
// Parameters:
// atom  - vector of atoms
// bond - vector of bonds
// n_bond - number of bonds
// dist - threshold bond length
void collapse_bonds(vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, double dist);

// Function: fix_one_sided_bonds()
//
// T-style bonds are corrected by creating an atom at the intersection point of T
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// thickness - bond thickness
// avg - average bond length
//
// Returns:
// new value for n_bond
int fix_one_sided_bonds(vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom, double thickness, double avg);

// Function: find_wedge_bonds()
//
// Detects wedge stereo-bonds
//
// Parameters:
// image -  original image
// atom - vector of atoms
// n_atom - number of atoms
// bond - vector of bonds
// n_bond - number of bonds
// bgColor - background color
// THRESHOLD_BOND - black-white binarization threshold
// max_dist_double_bond - maximum distance between double bonds
// avg - average bond length
// limit - minimum difference between the thick and thin ends of the bond
// dist - step away from the bond end to measure thickness
//
// Returns:
// Bond thickness
double find_wedge_bonds(const Image &image, vector<atom_t> &atom, int n_atom, vector<bond_t> &bond, int n_bond,
                        const ColorGray &bgColor, double THRESHOLD_BOND, double max_dist_double_bond, double avg, int limit, int dist = 0);

// Function: collapse_double_bonds()
//
// Removes a small bond connecting the ends of a conjoined double bond
//
// Parameters:
// bond  - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// dist - maximum bond length for small connecting bond
void collapse_double_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double dist);

// Function: find_up_down_bonds()
//
// Detects E/Z stereochemistry around a double bond
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom  - vector of atoms
// thickness - bond thickness
void find_up_down_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double thickness);

// Function: find_old_aromatic_bonds()
//
// Detects aromatic bonds which are drawn as a ring with a circle inside
//
// p - Potrace vectorization output
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// n_atom - number of atoms
// avg - average bond length
void find_old_aromatic_bonds(const potrace_path_t *p, vector<bond_t> &bond, int n_bond, vector<atom_t> &atom,int n_atom, double avg);

// Function: flatten_bonds()
//
// "Flattens" (combines into one bond) bonds which have a small "kink" in the middle
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// maxh - maximum vertical distance for the "kink"
void flatten_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double maxh);

// Function: mark_terminal_atoms()
//
// Finds and flags end atoms on a terminal bond
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom  - vector of atoms
// n_atom -  number of atoms
void mark_terminal_atoms(const vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, int n_atom);

// Function: find_limits_on_avg_bond()
//
// Calculates minimum and maximum bond length consistent with high confidence values
//
// Parameters:
// min_bond - a reference to minimum bond length
// max_bond - a reference to maximum bond_length
// pages_of_avg_bonds - a vector (element per page) of vectors (element per structure) of average bond lengths
// pages_of_ind_conf - a vector (pages) of vectors (structures on a page) of confidence estimes
void find_limits_on_avg_bond(double &best_bond, const vector<vector<double> > &pages_of_avg_bonds,
                             const vector<vector<double> > &pages_of_ind_conf);


// Function: find_wavy_bonds()
//
// Finds undefined stereo bonds (wavy bonds)
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// avg - average bond length
//
// Returns: number of bonds
int find_wavy_bonds(vector<bond_t> &bond,int n_bond,const vector<atom_t> &atom,double avg);


void remove_small_bonds_in_chars(vector<atom_t> &atom, vector<bond_t> &bond,vector<letters_t> &letters);
