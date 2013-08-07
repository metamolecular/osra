
#include <iostream>

#include <stdlib.h>
#include <string.h>

#include "ocradlib.h"

// Required by OCRAD, not used here:
#include <vector>
#include <stdio.h>

#include "common.h"
#include "rectangle.h"
#include "bitmap.h"
#include "blob.h"
#include "character.h"

using namespace std;

/* Actual max height is 12, but we leave some more for extensions: */
const char* TESTS[][12] =
{
  /* Test1: "N" is not detected */
  {
    "11100001",
    "11110001",
    "11110001",
    "11011001",
    "11011101",
    "11001111",
    "11000111",
    "11000111",
    "11000011",
  },
  /* Test2: "N" is not detected */
  {
    "11000011",
    "11100011",
    "11110011",
    "11110011",
    "10011011",
    "10011111",
    "10001111",
    "10000111",
    "10000111",
  },
  /* Test3: Detected as "r": */
  {
    "000000010",
    "000000111",
    "000001110",
    "000011100",
    "000111000",
    "001111000",
    "011110000",
    "111100000",
    "011000000",
  },
  /* Test4: Detected as "r": */
  {
    "00000011111",
    "00000111100",
    "00011110000",
    "01111100000",
    "11110000000",
    "11000000000",
  },
  /* Test5: Detected as "r": */
  {
    "00111",
    "01110",
    "01110",
    "01110",
    "11100",
    "11100",
    "11100",
    "11000",
    "11000",
  },
  /* Test6: Detected as "t": */
  {
    "00111",
    "00111",
    "00110",
    "01110",
    "01110",
    "01110",
    "11100",
    "11100",
    "11100"
  }
};

char run_test(int n)
{
  int height = 0;
  int width = strlen(TESTS[n][0]);

  while (TESTS[n][height] != NULL)
    {
      height++;
    }

  const char** image = TESTS[n];

  cout << "Test " << n + 1 << ": width x height = " << width << "x" << height << endl;

  // Blob for recognition attempt via Character::recognize1():
  Blob* b = new Blob(0, 0, width-1, height-1);

  // OCRAD_Pixmap for recognition attempt via OCRAD_result_first_character():
  struct OCRAD_Pixmap* opix = new OCRAD_Pixmap();

  unsigned char* bitmap_data = (unsigned char*) malloc(width * height);
  unsigned char* greymap_data = (unsigned char*) malloc(width * height);

  memset(bitmap_data, 0, width * height);
  memset(greymap_data, 255, width * height);

  opix->height = height;
  opix->width = width;

  opix->mode = OCRAD_bitmap;
  opix->data = bitmap_data;
//	opix->mode = OCRAD_greymap; opix->data = greymap_data;

  for (int row = 0; row < height; row++)
    {
      for (int col = 0; col < width; col++)
        if (image[row][col] == '1')
          {
            b->set_bit(row, col, true);
            bitmap_data[row * width + col] = 1;
            greymap_data[row * width + col] = 0;
          }
    }

  b->find_holes();

  Control control;
  Character a(b);

  // The Blob object was delegated to Character, which will free it on destruction:
  b = NULL;

  a.recognize1(control.charset, Rectangle::Rectangle(a.left(), a.top(), a.right(), a.bottom()));
  char c1 = a.byte_result();

  // Was the character recognised by OCRAD?
  cout << "+ recognised via Character::recognize1(): " << c1 << endl;

  char c2 = 0;
  OCRAD_Descriptor * const ocrdes = OCRAD_open();

  if (ocrdes && OCRAD_get_errno(ocrdes) == OCRAD_ok &&
      OCRAD_set_image(ocrdes, opix, 0) == 0 &&
      ( height >= 10 || OCRAD_scale( ocrdes, 2 ) == 0 ) &&
      OCRAD_recognize(ocrdes, 0) == 0 )
    c2 = OCRAD_result_first_character(ocrdes);

  OCRAD_close(ocrdes);

  delete opix;
  free(bitmap_data);
  free(greymap_data);

  cout << "+ recognised via OCRAD_result_first_character(): " << c2 << endl;
}

int main()
{
  for (unsigned int n = 0; n < 6; n++)
    {
      run_test(n);
    }
}
