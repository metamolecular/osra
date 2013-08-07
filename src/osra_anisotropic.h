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

// Header: osra_anisotropic.h
//
// Defines types and functions for anisotropic smoothing module.
//

#include <Magick++.h> // Magick::Image

//
// Section: Functions
//

// Function: anisotropic_smoothing()
//
// Performs Greycstoration anisotropic smoothing on an image according to the specified parameters
//
// Parameters:
//      image - image object
//      width - width of image
//      height - height of image
//      amplitude - amplitude of smoothing
//      sharpness - sharpness parameter
//      anisotropy - anisotropy parameter
//      alpha - alpha parameter for smoothing
//      sigma - sigma parameter for smoothing'
//
// Returns:
//      image object
//
// See also:
//      <anisotropic_scaling()>
Magick::Image anisotropic_smoothing(const Magick::Image &image, int width, int height, const float amplitude,
                                    const float sharpness, const float anisotropy, const float alpha, const float sigma);

// Function: anisotropic_scaling()
//
// Performs Greycstoration anisotropic scaling on an image
//
// Parameters:
//      image - image object
//      width - width of image
//      height - height of image
//      nw - new width
//      nh - new height
//
// Returns:
//      image object
//
// See also:
//      <anisotropic_smoothing()>
Magick::Image anisotropic_scaling(const Magick::Image &image, int width, int height, int nw, int nh);

