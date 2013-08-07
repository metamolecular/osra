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

#include <stddef.h> // NULL
#include <stdlib.h> // malloc(), free()

#include <iostream> // std::cout
#include <fstream> // std::ifstream
#include <sstream> // std:ostringstream

#include <osra_lib.h>

using namespace std;

int main(int argc, char **argv)
{
  if (argc < 2)
    {
      cout << "Usage: " << argv[0] << " [image_file_name]" << endl;
      return 1;
    }

  ifstream is(argv[1]);

  if (!is.is_open())
    {
      cout << "Failed to open a file '" << argv[1] << '\'' << endl;
      return 2;
    }

  // Learn the file size:
  is.seekg(0, ios::end);
  const int buf_size = (int) is.tellg();
  is.seekg(0, ios::beg);

  // Allocate memory:
  char* buf = (char*) malloc(buf_size);

  if (buf == NULL)
    {
      cout << "Failed to allocate " << buf_size << " bytes of memory" << endl;
      is.close();
      return 3;
    }

  is.read(buf, buf_size);

  // Call OSRA:
  const int result = osra_process_image(
                       buf,
                       buf_size,
                       cout,
                       0,
                       false,
                       0,
                       0,
                       0,
                       false,
                       false,
                       "sdf",
                       "",
                       true,
                       false,
                       false,
                       true,
                       true,
                       "",
                       "",
                       "",
                       false,
                       false
                     );

  // Release the allocated resources:
  is.close();
  free(buf);

  return result;
}
