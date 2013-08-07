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

#include <Magick++.h> // Magick::Image

// Header: unpaper.h
//
// Defines types and functions for unpaper image adjustment module.
//

//
// Section: Functions
//

// Function: unpaper()
//
// Performs unpaper image adjustment based on http://unpaper.berlios.de/
//
// Parameters:
//      picture - image object
//
// Returns:
//      0 in case of success or non-zero error code otherwise
int unpaper(Magick::Image &picture, double &radians, int &unpaper_dx, int &unpaper_dy);
