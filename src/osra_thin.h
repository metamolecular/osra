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

// Header: osra_thin.h
//
// Image thinning routines and noise factor computation
//
#include <Magick++.h>

using namespace Magick;

//
// Section: Functions
//

// Function: noise_factor()
//
// computes attributes of line thickness histogram
//
// Parameters:
// image - image to be processed
// width - image width
// height - image height
// bgColor - background color
// THRESHOLD_BOND - black-white binarization threshold
// resolution - resolution for which we're performing processing
// max - position of the maximum of the thickness histogram (most common thickness)
// nf45 - ratio of number of lines with thickness 4 to the number of lines with thickness 5
//
// Returns:
// Ratio of the number of lines with thickness 2 to number of lines of thickness 3
// or, if max == 2, ratio of the count of lines with thickness 1 to number of lines of thickness 2
// or, if max == 1, ratio of the count of lines with thickness 2 to number of lines of thickness 1
double noise_factor(const Image &image, int width, int height, const ColorGray &bgColor, double THRESHOLD_BOND,
                    int resolution, int &max, double &nf45);

// Function: thin_image()
//
// Performs image thinning based on Rosenfeld's algorithm
//
// Parameters:
// box - original image
// THRESHOLD_BOND - black-white binarization threshold
// bgColor - background color
//
// Returns:
// Thinned image
Image thin_image(const Image &box, double THRESHOLD_BOND, const ColorGray &bgColor);
