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

// Header: osra_labels.h
//
// declares functions dealing with atomic labels
//
#ifndef OSRA_LABELS_H
#define OSRA_LABELS_H

#include <string>
#include <vector>

#include <Magick++.h>

extern "C" {
#include <potracelib.h>
}

#include "osra.h"

using namespace std;

// struct: letters_s
// character found as part of atomic label
struct letters_s
{
  // double: x,y,r
  // coordinates of the center and radius of the circumscribing circle
  double x, y, r;
  // int: min_x, min_y, max_x, max_y
  // box coordinates
  int min_x, min_y,max_x,max_y;
  // char: a
  // character
  char a;
  // bool: free
  // whether or not it was already assign to an existing atomic label
  bool free;
  // pointer: curve
  //    pointer to the curve found by Potrace
  const potrace_path_t *curve;
};
// typedef: letters_t
// defines letters_t type based on letters_s struct
typedef struct letters_s letters_t;

// struct: label_s
// atomic label
struct label_s
{
  // doubles: x1,y1, x2, y2, r1, r2
  // central coordinates and circumradii for the first and last characters
  double x1, y1, r1, x2, y2, r2;
  // int: min_x, min_y, max_x, max_y
  // box coordinates
  int min_x, min_y,max_x,max_y;
  // string: a
  // atomic label string
  string a;
  // array: n
  // vector of character indices comprising the atomic label
  vector<int> n;
};
// typedef: label_t
// defines label_t type based on label_s struct
typedef struct label_s label_t;

//struct: lbond_s
//pairs of characters used for constucting atomic labels in <osra.cpp::assemble_labels()>
struct lbond_s
{
  //int: a,b
  // indices of first and second character in a pair
  int a, b;
  //double: x
  // x-coordinate of the first character
  double x;
  //bool: exists
  //pair of characters is available
  bool exists;
};
//typedef: lbond_t
//defines lbond_t type based on lbond_s struct
typedef struct lbond_s lbond_t;

//
// Section: Functions
//

// Function: assemble_labels()
//
// assembles characters into a string for a superatom label
//
// Parameters:
// letters - a vector of recognized characters
// n_letters - the number of recognized characters
// label - a reference to the vector which will contain superatom labels
//
// Returns:
// number of superatom labels
int assemble_labels(vector<letters_t> &letters, int n_letters, vector<label_t> &label);

// Function: find_chars()
//
// searches for perspective characters in the image and calls OCR routines
//
// Parameters:
// p - vectorized output of Potrace routines
// orig - original image
// letters - vector which will contain recognized characters
// atom - vector of atoms
// bond - vector of bonds
// n_atom - number of atoms
// n_bond - number of bonds
// height - image height
// width - image width
// bgColor - background color
// THRESHOLD - black-white binarization threshold
// max_font_width - maximum font width for the specific resolution in pixels
// max_font_height - maximum font height for the specific resolution in pixels
// real_font_width - detected font width
// real_font_height - detected font height
// verbose - flag for verbose output
//
// Returns:
// number of recognized characters
int find_chars(const potrace_path_t * p, const Image &orig, vector<letters_t> &letters, vector<atom_t> &atom, vector<
               bond_t> &bond, int n_atom, int n_bond, int height, int width, ColorGray &bgColor, double THRESHOLD,
               int max_font_width, int max_font_height, int &real_font_width, int &real_font_height, bool verbose);

// Function: find_numbers()
//
// searches for numbers 0..9 in the image and calls OCR routines
//
// Parameters:
// p - vectorized output of Potrace routines
// orig - original image
// letters - vector which will contain recognized characters
// atom - vector of atoms
// bond - vector of bonds
// n_atom - number of atoms
// n_bond - number of bonds
// height - image height
// width - image width
// bgColor - background color
// THRESHOLD - black-white binarization threshold
// n_letters - number of previously recognized characters
//
// Returns:
// number of recognized characters
int find_numbers(const potrace_path_t * p, const Image &orig, vector<letters_t> &letters, vector<atom_t> &atom, vector<bond_t> &bond, 
		 int n_atom, int n_bond, int height, int width, ColorGray &bgColor, double THRESHOLD, int n_letters);

// Function: find_plus_minus()
//
// Detects plus and minus signs in atomic charge labels
//
// Parameters:
// p - Potrace vectorization output
// image - original image
// bgColor - background color
// THRESHOLD - black-white binarization threshold
// letters - a vector of atomic label characters
// atom -  a vector of atoms
// bond - a vector of bonds
// n_atom - number of atoms
// n_bond - number of bonds
// height - image height
// width - image width
// max_font_height - maximum font height for a given resolution
// max_font_width - maximum font width for a given resolution
// n_letters - number of characters
//
// Returns:
// new number of characters
int find_plus_minus(const potrace_path_t *p, const Image &image, ColorGray &bgColor, double THRESHOLD, vector<letters_t> &letters, vector<atom_t> &atom, vector<bond_t> &bond,
                    int n_atom, int n_bond, int height, int width, int max_font_height, int max_font_width, int n_letters, double avg_bond_length);

// Function: clean_unrecognized_characters()
//
// Attempts to assemble collections of small bonds which did not pass OCR routines
// into an unspecified atomic label "*"
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// real_font_height - font height as determined by preceeding OCR routines
// real_font_width - font width as determined by preceeding OCR routines
// size - minimum number of bonds to constitute a perspective character
// letters - vector of recognized characters from OCR
// n_letters - number of recognized characters
//
// Returns:
// new value for n_letters
int clean_unrecognized_characters(vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom, int real_font_height,
                                  int real_font_width, unsigned int size, vector<letters_t> &letters, int n_letters);

// Function: remove_small_terminal_bonds()
//
// Small terminal bonds are often a result of unrecognized atomic labels
// This function replaces such bonds with "Xx" atomic label
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// avg - average bond length
void remove_small_terminal_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double avg);

// Function: remove_small_bonds()
//
// Removes very small single bonds or replaces them, if stand-alone and next to a character and vertical
// with a character "l"
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// letters - vector of characters
// n_letters  - number of characters
// max_font_height - maximum font height
// min_font_height - minimum font height
// avg - average bond length
//
// Returns:
// new value for n_letters
int remove_small_bonds(vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom, vector<letters_t> &letters,
                       int n_letters, int max_font_height, int min_font_height, double avg);

// Function: find_fused_chars()
//
// Attempts to recognize characters fused to a bond
//
// Parameters:
// bond - vector of bonds
// n_bond - number of bonds
// atom - vector of atoms
// letters - vector of characters
// n_letters - number of characters
// max_font_height - maximum font height
// max_font_width - maximum font width
// dummy - if set, character to substitute for OCR results
// orig  - original image
// bgColor - background color
// THRESHOLD - black-white threshold for image binarization
// size - minimum number of bonds which can constitute a character
// verbose - flag for verbose output
//
// Returns:
// new value for n_letters
int find_fused_chars(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, vector<letters_t> &letters, int n_letters,
                     int max_font_height, int max_font_width, char dummy, const Image &orig, const ColorGray &bgColor,
                     double THRESHOLD, unsigned int size, bool verbose);
#endif
