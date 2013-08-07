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

#include <stddef.h> // NULL
#include <stdlib.h> // free()
#include <string.h> // strlen()

#include <iostream> // std::cout

#include <tesseract/baseapi.h>

using namespace std;

const char* PICTURE[][50] =
{
  {
    "##.........",
    "##.........",
    "...........",
    "...........",
    "..........#",
    "..........#",
    "........###",
    "......#####",
    "....##.####",
    "...########",
    "###########"
  },
  {
    "##.........",
    "##.........",
    "...........",
    "...........",
    "..........#",
    "..........#",
    "........###",
    "......#####",
    "....##.####",
    "...########",
    "###########"
  },
  {
    "##.........",
    ".#.........",
    "...........",
    "...........",
    "..........#",
    "..........#",
    "..........#",
    "......#####",
    ".....######",
    "....#######",
    "###########"
  },
  {
    "...........",
    "...........",
    "...........",
    "...........",
    "..........#",
    "........#.#",
    "........###",
    "......#####",
    "..#########",
    ".##########",
    ".##########"
  }
};
// Create a binary pixel image and recognize it using the given engine:
void recognize(int n, tesseract::TessBaseAPI &t)
{
  int height = 0;
  int width = strlen(PICTURE[n][0]);

  while (PICTURE[n][height] != NULL)
    height++;

  cout << "Picture " << n + 1 << ": width x height = " << width << "x" << height << endl;

  unsigned char *pixmap = (unsigned char *) malloc(width * height);

  for (int row = 0; row < height; row++)
    for (int col = 0; col < width; col++)
      pixmap[row * width + col] = PICTURE[n][row][col] == '#' ? 255 : 0;

  char* text = t.TesseractRect(pixmap, 1, width, 0, 0, width, height);

  cout << "Result: " << text << '.' << endl;

  free(text);
  free(pixmap);
}

#ifdef TESS_GLOBAL_INSTANCE
// Global variable:
tesseract::TessBaseAPI tess;

int main(int argc, char **argv)
{
  tess.Init(NULL, "eng", NULL, 0, false);

  for (unsigned int n = 0; n < sizeof(PICTURE) / 200; n++)
    {
      recognize(n, tess);

      tess.Clear();
    }

  tess.End();
}
#else
int main(int argc, char **argv)
{
  for (unsigned int n = 0; n < sizeof(PICTURE) / 200; n++)
    {
      tesseract::TessBaseAPI tess;
      tess.Init(NULL, "eng", NULL, 0, false);

      recognize(n, tess);

      tess.End();
    }
}
#endif
