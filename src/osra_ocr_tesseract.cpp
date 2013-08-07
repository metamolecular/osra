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
#include <stdlib.h> // free()
#include <ctype.h> // isalnum()
#include <string.h> // strlen()

#include <string> // std::string

#include <tesseract/baseapi.h>

const char UNKNOWN_CHAR = '_';

using namespace std;

// Global variable:
tesseract::TessBaseAPI tess;

void osra_tesseract_init()
{
  tess.Init(NULL, "eng", tesseract::OEM_DEFAULT, NULL, 0, NULL, NULL, false);
}

void osra_tesseract_destroy()
{
  tess.End();
}

char osra_tesseract_ocr(unsigned char *pixmap, int width, int height, const string &char_filter)
{
  char result = UNKNOWN_CHAR;

  char *text = tess.TesseractRect(pixmap, 1, width, 0, 0, width, height);

  // TODO: Why text length should be exactly 3? Give examples...
  if (text != NULL && strlen(text) == 3 && isalnum(text[0]) && (char_filter.empty() || char_filter.find(text[0], 0) != string::npos))
    result = text[0];

  free(text);

  return result;
}
