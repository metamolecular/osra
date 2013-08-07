/******************************************************************************
 OSRA: Optical Structure Recognition

 This is a U.S. Government work (2007-2010) and is therefore not subject to
 copyright. However, portions of this work were obtained from a GPL or
 GPL-compatible source.
 Created by Igor Filippov, 2007-2010 (igorf@helix.nih.gov)

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

#include <Magick++.h>

#include <string>
#include <iostream>

using namespace std;

int main(int argc, char **argv)
{
  string fileName = argc > 1 ? argv[1] : "aaaa";
  string type;

  Magick::InitializeMagick(*argv);

  try
    {
      Magick::Image image;
      image.ping(fileName);
      type = image.magick();
    }
  catch (...)
    {
      cerr << "Cannot open file '" << fileName << "'" << endl;
      exit(1);
    }

  cerr << "File type is '" << type << "'" << endl;

  return 0;
}
