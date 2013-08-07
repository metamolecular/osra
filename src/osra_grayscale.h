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

// Header: osra_grayscale.h
//
// Declares grayscale conversion functions
//

#include <Magick++.h>

using namespace std;
using namespace Magick;

//
// Section: Functions
//

// Function: getBgColor()
//
// Detects the backgroun color of the image
//
// Parameters:
// image -  a reference to the image object
//
// Returns:
// a Color object corresponding to the background color
const Color getBgColor(const Image &image);

// Function: convert_to_gray()
//
// Converts image to grayscale
//
// Parameters:
// image - reference to Image object
// invert - flag set if the image is white-on-black
// adaptive - flag set if adaptive thresholding is enforced
// verbose - flag set if verbose reporting is on
//
// Returns:
// a boolean flag indicating whether adaptive thresholding is indicated
bool convert_to_gray(Image &image, bool invert, bool adaptive, bool verbose);
