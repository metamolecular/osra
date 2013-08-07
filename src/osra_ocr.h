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

// Header: osra_ocr.h
//
// Defines types and functions for OSRA OCR module.
//

#include <string> // std::string
#include <map> // std::map

#include <Magick++.h> // Magick::Image, Magick::ColorGray

using namespace std;

//
// Section: Functions
//

// Function: osra_ocr_init()
//
// Initialises OCR engine. Should be called at e.g. program startup.
//
void osra_ocr_init();

// Function: osra_ocr_destroy()
//
// Releases all resources allocated by OCR engine.
//
void osra_ocr_destroy();

// Function: get_atom_label()
//
// OCR engine function, does single character recognition
//
// Parameters:
//      image - image object
//      bg - gray-level background color
//      x1, y1, x2, y2 - coordinates of the character box
//      THRESHOLD - graylevel threshold for image binarization
//      dropx, dropy - coordinates of drop point from where breadth-first algorithm will search for single connected component
//                     which is hopefully the character we are trying to recognize
//      no_filtering - do not apply character filter
//      numbers - only allow numbers in the output 0..9
//
// Returns:
//      recognized character or 0
char get_atom_label(const Magick::Image &image, const Magick::ColorGray &bg, int x1, int y1, int x2, int y2,
                    double THRESHOLD, int dropx, int dropy, bool no_filtering, bool verbose, bool numbers = false);

// Function: fix_atom_name()
//
// Corrects common OCR errors by using spelling dictionary
//
// Parameters:
//      s - Original atomic label as returned by OCR engine.
//      n - The number of bonds attached to the atom.
//      fix - spelling dictionary
//      superatom - dictionary of superatom labels mapped to SMILES
//      debug - enables output of debugging information to stdout
//
// Returns:
//      corrected atomic label
const string fix_atom_name(const string &s, int n, const map<string, string> &fix,
                           const map<string, string> &superatom, bool debug);
