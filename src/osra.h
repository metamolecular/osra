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

// Header: osra.h
//
// Defines types and functions exported from main module to other modules.
//
#ifndef OSRA_H
#define OSRA_H

#include <string> // std:string
#include <vector> // std::vector

#include <Magick++.h> // Magick::Image, Magick::ColorGray

extern "C" {
#include <potracelib.h>
}

using namespace std;
using namespace Magick;

// struct: atom_s
//      Contains information about perspective atom
struct atom_s
{
atom_s(double xx=0, double yy=0, const potrace_path_t* p=NULL) : 
  x(xx),y(yy),min_x(xx),min_y(yy),max_x(xx),max_y(yy),curve(p),label(" "),n(0),anum(0), exists(false),corner(false),terminal(false),charge(0) {}
  // doubles: x, y
  //    coordinates within the image clip
  double x, y;
  // string: label
  //    atomic label
  string label;
  // int: n
  //    counter of created OBAtom objects in <create_molecule()>
  int n;
  // int: anum
  //    atomic number
  int anum;
  // pointer: curve
  //    pointer to the curve found by Potrace
  const potrace_path_t *curve;
  // bools: exists, corner, terminal
  //    atom exists, atom is at the corner (has two bonds leading to it), atom is a terminal atom
  bool exists, corner, terminal;
  // int: charge
  //    electric charge on the atom
  int charge;
  // int: min_x, min_y, max_x, max_y
  // box coordinates
  int min_x, min_y,max_x,max_y;
};
// typedef: atom_t
//      defines atom_t type based on atom_s struct
typedef struct atom_s atom_t;

// struct: bond_s
//      contains information about perspective bond between two atoms
struct bond_s
{
bond_s(int i=0, int j=0, const potrace_path_t* p=NULL) : 
  a(i),b(j),curve(p),type(1),exists(true),hash(false),wedge(false),up(false),down(false),Small(false),arom(false),conjoined(false) {}
  // ints: a, b, type
  //    starting atom, ending atom, bond type (1=single, 2=double, 3=triple)
  int a, b, type;
  // pointer: curve
  //    pointer to the curve found by Potrace
  const potrace_path_t *curve;
  // bools: exists, hash, wedge, up, down, Small, arom
  //    bond existence and type flags
  bool exists;
  bool hash;
  bool wedge;
  bool up;
  bool down;
  bool Small;
  bool arom;
  // bool: conjoined
  //    true for a double bond which is joined at one end on the image
  bool conjoined;
};
// typedef: bond_t
//      defines bond_t type based on bond_s struct
typedef struct bond_s bond_t;

// Section: Constants
//
// Constants: global defines
//
// MAX_ATOMS - maximum size of the vector holding perspective atoms
// MAX_FONT_HEIGHT - maximum font height at a resolution of 150 dpi
// MAX_FONT_WIDTH - maximum font width at a resolution of 150 dpi
// MIN_FONT_HEIGHT - minimum font height
// BG_PICK_POINTS - number of points to randomly pick to determine background color
// D_T_TOLERANCE - cosine tolerance to find parallel bonds for double-triple bond extraction
// V_DISPLACEMENT - threshold vertical displacement in pixels
// DIR_CHANGE - threshold direction change in pixels
// THRESHOLD_GLOBAL - gray-level threshold for image binarization
// THRESHOLD_LOW_RES - gray-level threshold for low resolutions (72 dpi)
// MAX_RATIO - maximum black/white fill ratio for perspective molecular structures
// MIN_ASPECT - minimum aspect ration
// MAX_ASPECT - maximum aspect ratio
// MIN_A_COUNT - minimum number of atoms
// MAX_A_COUNT - maximum number of atoms
// MIN_CHAR_POINTS - minimum number of black and white pixels in a character box
// MAX_BOND_THICKNESS - maximum bond thickness
// SMALL_PICTURE_AREA - threshold area of the image to be consider a small picture
// NUM_RESOLUTIONS - number of resolutions to try
// MAX_DASH - maximum size of a dash in a dashed bond
// CC_BOND_LENGTH - average carbon-carbon bond length
// FRAME - border around structure in a segmented image
// SEPARATOR_ASPECT - aspect ratio for a perspective separator line
// SEPARATOR_AREA - area for a perspective separator line
// MAX_DIST - maximum distance in pixels between neighboring segments in image segmentation routines
// MAX_AREA_RATIO - maximum area ratio for connected compoments in image segmentation
// SINGLE_IMAGE_DIST - default distance between connected components in a single structure image
// THRESHOLD_LEVEL - threshold level for feature matrix for image segmentation
// TEXT_LINE_SIZE - maximum atomic label size in characters
// PARTS_IN_MARGIN - take only every other pixel on a connected component margin for speed
// BORDER_COUNT - threshold number of pixels on a box border to be considered a table
// MAX_SEGMENTS - maximum number of connected compoment segments
// MAX_FRAGMENTS - maximum number of fragments
// STRUCTURE_COUNT - threshold number of structures to compute limits on average bond length
// SPELLING_TXT - spelling file for OCR corrections
// SUPERATOM_TXT - superatom file for mapping labels to SMILES
#define PI 3.14159265358979323846
#define MAX_ATOMS 10000
#define MAX_FONT_HEIGHT 22
#define MAX_FONT_WIDTH 21
#define MIN_FONT_HEIGHT 5
#define BG_PICK_POINTS 1000
#define D_T_TOLERANCE 0.95
#define V_DISPLACEMENT 3
#define DIR_CHANGE 2
#define THRESHOLD_GLOBAL 0.4
#define THRESHOLD_LOW_RES 0.2
#define MAX_RATIO 0.2
#define MIN_ASPECT 0.1
#define MAX_ASPECT 10.
#define MIN_A_COUNT 5
#define MAX_A_COUNT 250
#define MIN_CHAR_POINTS 2
#define MAX_BOND_THICKNESS 10
#define SMALL_PICTURE_AREA 6000
#define NUM_RESOLUTIONS 5
#define MAX_DASH 80
#define CC_BOND_LENGTH 1.5120
#define FRAME 5
#define SEPARATOR_ASPECT 100
#define SEPARATOR_AREA 300
#define MAX_DIST 50
#define MAX_AREA_RATIO 50
#define SINGLE_IMAGE_DIST 1000
#define MAX_DISTANCE_BETWEEN_ARROWS 200
#define THRESHOLD_LEVEL 4
#define TEXT_LINE_SIZE 8
#define PARTS_IN_MARGIN 2
#define BORDER_COUNT 3000
#define MAX_SEGMENTS 10000
#define MAX_FRAGMENTS 10
#define STRUCTURE_COUNT 20
#define SPELLING_TXT "spelling.txt"
#define SUPERATOM_TXT "superatom.txt"
#define RECOGNIZED_CHARS "oOcCnNHFsSBuUgMeEXYZRPp23456789AmThD"

#define ERROR_SPELLING_FILE_IS_MISSING          -1
#define ERROR_SUPERATOM_FILE_IS_MISSING         -2
#define ERROR_OUTPUT_FILE_OPEN_FAILED           -3
// This error code may be returned, if ImageMagic was not able to find the .mgk files.
// Check that MAGICK_CONFIGURE_PATH points to the location of *.mgk configuration files (check here http://www.imagemagick.org/script/resources.php).
#define ERROR_UNKNOWN_IMAGE_TYPE                -4
#define ERROR_ILLEGAL_ARGUMENT_COMBINATION      -5
// This error code usually means:
// (a) You have no /usr/lib/openbabel/x.x.x/smilesformat.so library installed. Install the format libraries / check http://openbabel.org/docs/dev/Installation/install.html#environment-variables
// (b) The format libraries are installed, but do not correspond to /usr/lib/libopenbabel.so.y.y.y. Check they correspond to the same OpenBabel version.
// (c) You need to preload OpenBabel e.g. using LD_PRELOAD=/usr/lib/libopenbabel.so
#define ERROR_UNKNOWN_OPENBABEL_FORMAT          -6

#endif
