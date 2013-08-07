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

// Header: osra_common.h
//
// Common functionality routines
//

extern "C" {
#include <potracelib.h>
}

#include "osra.h"
#include "osra_segment.h"

#include <vector> // std::vector

// Defines necessary for potrace vectorization functions
#define BM_WORDSIZE ((int)sizeof(potrace_word))
#define BM_WORDBITS (8*BM_WORDSIZE)
#define BM_HIBIT (((potrace_word)1)<<(BM_WORDBITS-1))
#define bm_scanline(bm, y) ((bm)->map + (y)*(bm)->dy)
#define bm_index(bm, x, y) (&bm_scanline(bm, y)[(x)/BM_WORDBITS])
#define bm_mask(x) (BM_HIBIT >> ((x) & (BM_WORDBITS-1)))
#define bm_range(x, a) ((int)(x) >= 0 && (int)(x) < (a))
#define bm_safe(bm, x, y) (bm_range(x, (bm)->w) && bm_range(y, (bm)->h))
#define BM_USET(bm, x, y) (*bm_index(bm, x, y) |= bm_mask(x))
#define BM_UCLR(bm, x, y) (*bm_index(bm, x, y) &= ~bm_mask(x))
#define BM_UPUT(bm, x, y, b) ((b) ? BM_USET(bm, x, y) : BM_UCLR(bm, x, y))
#define BM_PUT(bm, x, y, b) (bm_safe(bm, x, y) ? BM_UPUT(bm, x, y, b) : 0)

//
// Section: Functions
//

// Function: get_pixel()
//
// Returns a binarized pixel value from a gray-level image
//
// Parameters:
//      image - image object
//      bg - gray-level background color
//      x, y - coordinates of the pixel
//      THRESHOLD - gray-level threshold for binarization
//
// Returns:
//      1 for set pixel, 0 for background
int get_pixel(const Magick::Image &image, const Magick::ColorGray &bg, unsigned int x, unsigned int y, double THRESHOLD);

// Function: trim()
//
// Remove leading and trailing whitespace
//
// Parameters:
//      s - string to trim (in/out parameter)
void trim(std::string &s);

// Function: distance()
//
// Returns L2 measure between 2 points in 2d  space
//
// Parameters:
// x1, y1 - coordinates of the first point
// x2,y2 - coordinates of the second point
//
// Returns:
// Distance between points
double distance(double x1, double y1, double x2, double y2);

// Function: bond_length()
//
// Returns the length of i-th bond
//
// Parameters:
// bond - vector of bond_t objects in molecule
// i - bond number
// atom - vector of atom_t objects in molecule
//
// Returns:
// Length of i-th bond
double bond_length(const vector<bond_t> &bond, int i, const vector<atom_t> &atom);

// Function: delete_curve()
//
// Deletes atoms and bonds which belong to a give curve
//
// Parameters:
// atom - vector of atom_t objects in molecule
// bond - vector of bond_t objects in molecule
// n_atom - last atom in molecule plus one
// n_bond last bond in molecule plus one
// curve - potrace_path_t pointer denoting a curve which should be deleted
//
// Returns:
// nothing
void delete_curve(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond,
                  const potrace_path_t * const curve);

// Function: delete_curve_with_children()
//
// Same as delete_curve() but also removes all the inner curves
//
// Parameters:
// atom - vector of atom_t objects in molecule
// bond - vector of bond_t objects in molecule
// n_atom - last atom in molecule plus one
// n_bond last bond in molecule plus one
// curve - potrace_path_t pointer denoting a curve which should be deleted
//
// Returns:
// nothing
void delete_curve_with_children(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond,
                                const potrace_path_t * const p);

// Function: detect_curve()
//
// Detects if a given curve correspond to any single, simple, non-deleted bonds
//
// Parameters:
// bond - vector of bonds in molecule
// n_bond - number of all bonds (last bond plus one)
// curve - curve pointer we are checking
//
// Returns:
// True if a single bond corresponding to such curve is found, False otherwise
bool detect_curve(vector<bond_t> &bond, int n_bond, const potrace_path_t * const curve);

// Function: terminal_bond()
//
// Checks if there exists a bond which has a specified atom as one of the ends and is not a given bond
//
// Parameters:
// a - atom number to check if it belongs to a bond
// b - bond number to avoid
// bond - vector of bonds in a molecule
// n_bond - number of bonds in a molecule
//
// Returns:
// False if a bond ending on atom a is found, True otherwise
bool terminal_bond(int a, int b, const vector<bond_t> &bond, int n_bond);

// Function: bm_new()
//
// Reserves space for Potrace object for image vectorization
//
// Parameters:
// w - width, h - height
//
// Returns:
// pointer to potrace_bitmap_t
potrace_bitmap_t *const bm_new(int w, int h);

// Function: angle4()
//
// Returns cosine of the angle between two vectors given by their end points
//
// Parameters:
// x1,y1, x2,y2 - coordinates of the end points of the first vector
// x3,y3, x4, y4 - coordinates of the end points of the second vector
//
// Returns:
// cosine of the angle between two vectors
double angle4(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

// Function: angle_between_bonds()
//
// Returns cosine of the angle between i-th and j-th bonds in a molecule
//
// Parameters:
// bond - vector of bonds in a molecule
// i, j - bonds for which we are calculate the angle
// atom - vector of atoms in a molecule
//
// Returns:
// cosine of the angle between bonds i and j
double angle_between_bonds(const vector<bond_t> &bond, int i, int j, const vector<atom_t> &atom);

// Function: distance_from_bond_y()
//
// Returns the length of a normal vector dropped from a point to a line containing the bond
//
// Parameters:
// x0,y0,x1,y1 - coordinates of the bond ends
// x,y - coordinates of the point for which we are looking for the distance from the bond
//
// Returns:
// Vertical distance from a point to the bond axis if the coordinates are rotated so that the bond is horizontal
double distance_from_bond_y(double x0, double y0, double x1, double y1, double x, double y);

// Function: distance_between_bonds()
//
// Calculates maximum vertical distance from j-th bond ends to bond i.
//
// Parameters:
// bond - vector of bonds in a molecule
// i, j - bonds
// atom - vector of atoms in a molecule
//
// Return:
// Maximum vertical distance from bond j ends to bond i when the coordinates are rotated so that bond i
// is horizontal. If the difference between distances from end a and end b is greater than 4 FLT_MAX
// is returned.
double distance_between_bonds(const vector<bond_t> &bond, int i, int j, const vector<atom_t> &atom);

// Function: distance_from_bond_x_a()
//
// Horizontal distance from end a of a bond to a specified point
//
// Parameters:
// x0,y0,x1,y1 - coordinates of bond ends
// x, y - coordinates of the specified point
//
// Returns:
// Horizontal component of the distance between point x,y and end a of a bond
// when coordinates are rotated so that the bond is horizontal. a->b is the positive direction.
double distance_from_bond_x_a(double x0, double y0, double x1, double y1, double x, double y);

// Function: distance_from_bond_x_b()
//
// Horizontal distance from end b of a bond to a specified point
//
// Parameters:
// x0,y0,x1,y1 - coordinates of bond ends
// x, y - coordinates of the specified point
//
// Returns:
// Horizontal component of the distance between point x,y and end b of a bond
// when coordinates are rotated so that the bond is horizontal. a->b is the positive direction
double distance_from_bond_x_b(double x0, double y0, double x1, double y1, double x, double y);

// Function: percentile75()
//
// Returns the bond length at 75% (third quntile)
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
//
// Returns:
// Third quintile value for bond length.
double percentile75(const vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom);

// Function: count_pages()
//
// Counts the number of pages in a specified image file
//
// Parameters:
// input - File name
//
// Returns:
// number of pages in an image file
int count_pages(const string &input);

// Function: count_atoms()
//
// Counts the number of non-deleted atoms in a molecule
//
// Parameters:
// atom - vector of atoms
// n_atom - number of atoms
//
// Returns:
// number of non-deleted atoms
int count_atoms(const vector<atom_t> &atom, int n_atom);

// Function: count_bonds()
//
// Counts the number of non-deleted bonds in a molecule
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
//
// Returns:
// number of non-deleted bonds
int count_bonds(const vector<bond_t> &bond, int n_bond, int &bond_max_type);

// Function: load_config_map()
//
// Loads text file into a std:map structure. Used for loading superatom and spelling files.
//
// Parameters:
// file - file name
// out - std:map object with the result
//
// Returns:
// True if file load is successful, False otherwise
bool load_config_map(const string &file, map<string, string> &out);

// Function: comp_boxes()
//
// Box coordinate comparison function used for sorting molecule-containing boxes
// in a top-down left-to-right order
//
// Parameters:
// aa,bb - box_t objects to compare
//
// Returns:
// True if box aa is higher or to the left of box bb.
// False otherwise
bool comp_boxes(const box_t &aa, const box_t &bb);

void debug_image(Image image, const vector<atom_t> &atom, int n_atom, const vector<bond_t> &bond, int n_bond,
		 const string &fname);
