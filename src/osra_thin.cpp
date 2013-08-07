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

#include "osra_common.h"
#include "osra_thin.h"
#include <iostream>

/*------------------- ThinImage - Thin binary image. --------------------------- *
 *
 *    Description:
 *        Thins the supplied binary image using Rosenfeld's parallel
 *        thinning algorithm.
 *
 *    On Entry:
 *        image = Image to thin.
 *
 *------------------------------------------------------------------------------- */

/* Direction masks:                  */
/*   N     S     W        E            */
static unsigned int masks[] = { 0200, 0002, 0040, 0010 };

/*    True if pixel neighbor map indicates the pixel is 8-simple and  */
/*    not an end point and thus can be deleted.  The neighborhood     */
/*    map is defined as an integer of bits abcdefghi with a non-zero  */
/*    bit representing a non-zero pixel.  The bit assignment for the  */
/*    neighborhood is:                                                */
/*                                                                    */
/*                            a b c                                   */
/*                            d e f                                   */
/*                            g h i                                   */

static unsigned char todelete[512] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1,
                                       1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                       0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1,
                                       1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                                       1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                       0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
                                       1, 1, 1, 1
                                     };

void thin1(unsigned char * const ptr, unsigned int xsize, unsigned int ysize)
{
  unsigned char *y_ptr, *y1_ptr;
  unsigned char bg_color = 0, colour = 1;
  unsigned int x, y; /* Pixel location               */
  unsigned int i; /* Pass index           */
  unsigned int pc = 0; /* Pass count           */
  unsigned int count = 1; /* Deleted pixel count          */
  unsigned int p, q; /* Neighborhood maps of adjacent*/
  /* cells                        */
  unsigned char *qb; /* Neighborhood maps of previous*/
  /* scanline                     */
  unsigned int m; /* Deletion direction mask      */

  qb = (unsigned char*) malloc(xsize * sizeof(unsigned char));
  qb[xsize - 1] = 0; /* Used for lower-right pixel   */

  while (count)   /* Scan image while deletions   */
    {
      pc++;
      count = 0;

      for (i = 0; i < 4; i++)
        {

          m = masks[i];

          /* Build initial previous scan buffer.                  */
          p = (ptr[0] == colour);
          for (x = 0; x < xsize - 1; x++)
            qb[x] = (unsigned char) (p = ((p << 1) & 0006) | (unsigned int) (ptr[x + 1] == colour));

          /* Scan image for pixel deletion candidates.            */
          y_ptr = ptr;
          y1_ptr = ptr + xsize;
          for (y = 0; y < ysize - 1; y++, y_ptr += xsize, y1_ptr += xsize)
            {
              q = qb[0];
              p = ((q << 2) & 0330) | (y1_ptr[0] == colour);

              for (x = 0; x < xsize - 1; x++)
                {
                  q = qb[x];
                  p = ((p << 1) & 0666) | ((q << 3) & 0110) | (unsigned int) (y1_ptr[x + 1] == colour);
                  qb[x] = (unsigned char) p;
                  if (((p & m) == 0) && todelete[p])
                    {
                      count++;
                      y_ptr[x] = bg_color; /* delete the pixel */
                    }
                }

              /* Process right edge pixel.                        */
              p = (p << 1) & 0666;
              if ((p & m) == 0 && todelete[p])
                {
                  count++;
                  y_ptr[xsize - 1] = bg_color;
                }
            }

          /* Process bottom scan line.                            */
          q = qb[0];
          p = ((q << 2) & 0330);

          y_ptr = ptr + xsize * (ysize - 1);
          for (x = 0; x < xsize; x++)
            {
              q = qb[x];
              p = ((p << 1) & 0666) | ((q << 3) & 0110);
              if ((p & m) == 0 && todelete[p])
                {
                  count++;
                  y_ptr[x] = bg_color;
                }
            }
        }
    }
  free(qb);
}

Image thin_image(const Image &box, double THRESHOLD_BOND, const ColorGray &bgColor)
{
  Image image(Geometry(box.columns(), box.rows()), "white");
  image.type(GrayscaleType);
  unsigned int xsize = box.columns();
  unsigned int ysize = box.rows();
  unsigned char *ptr = (unsigned char*) malloc(xsize * ysize * sizeof(unsigned char));

  for (unsigned int i = 0; i < xsize; i++)
    for (unsigned int j = 0; j < ysize; j++)
      ptr[i + j * xsize] = get_pixel(box, bgColor, i, j, THRESHOLD_BOND);

  if (xsize>1 && ysize>1)
    thin1(ptr, xsize, ysize);

  for (unsigned int i = 0; i < xsize; i++)
    for (unsigned int j = 0; j < ysize; j++)
      if (ptr[i + j * xsize] == 1)
        image.pixelColor(i, j, "black");

  free(ptr);
  return (image);
}


double noise_factor(const Image &image, int width, int height, const ColorGray &bgColor, double THRESHOLD_BOND,
                    int resolution, int &max, double &nf45)
{
  int max_thick = 40;
  vector<double> n(max_thick, 0);
  double nf;
  vector<int> lines;
  for (int i = 0; i < width; i++)
    {
      int j = 0;
      while (j < height)
        {
          while (!get_pixel(image, bgColor, i, j, THRESHOLD_BOND) && j < height)
            j++;
          int l = 0;
          while (get_pixel(image, bgColor, i, j, THRESHOLD_BOND) && j < height)
            {
              l++;
              j++;
            }
          if (l < max_thick)
            n[l]++;
	  lines.push_back(l);
        }
    }
  for (int i = 0; i < height; i++)
    {
      int j = 0;
      while (j < width)
        {
          while (!get_pixel(image, bgColor, j, i, THRESHOLD_BOND) && j < width)
            j++;
          int l = 0;
          while (get_pixel(image, bgColor, j, i, THRESHOLD_BOND) && j < width)
            {
              l++;
              j++;
            }
          if (l < max_thick)
            n[l]++;
	  lines.push_back(l);
        }
    }
  /*for (int i = 0; i < height; i++)
    {
      int j = 0;
      while (j < width && i+j < height)
        {
          while (!get_pixel(image, bgColor, j, i+j, THRESHOLD_BOND) && j < width && i+j < height)
            j++;
          int l = 0;
          while (get_pixel(image, bgColor, j, i+j, THRESHOLD_BOND) && j < width && i+j < height)
            {
              l++;
              j++;
            }
          if (l < max_thick)
            n[l]++;
	  lines.push_back(l);
        }
    }
  for (int i = 0; i < height; i++)
    {
      int j = width - 1;
      while (j >= 0 && i+(width-1-j) < height)
        {
          while (!get_pixel(image, bgColor, j, i+(width-1-j), THRESHOLD_BOND) && j >= 0 && i+(width-1-j) < height)
            j--;
          int l = 0;
          while (get_pixel(image, bgColor, j, i+(width-1-j), THRESHOLD_BOND) && j >= 0 && i+(width-1-j) < height)
            {
              l++;
              j--;
            }
          if (l < max_thick)
            n[l]++;
	  lines.push_back(l);
        }
    }
  for (int i = 0; i < width; i++)
    {
      int j = 0;
      while (j < height && i+j < width)
        {
          while (!get_pixel(image, bgColor, i+j, j, THRESHOLD_BOND) && j < height && i+j < width)
            j++;
          int l = 0;
          while (get_pixel(image, bgColor, i+j, j, THRESHOLD_BOND) && j < height && i+j < width)
            {
              l++;
              j++;
            }
          if (l < max_thick)
            n[l]++;
	  lines.push_back(l);
        }
    }
  for (int i = 0; i < width; i++)
    {
      int j = height - 1;
      while (j > 0 && i+(height-1-j) < width)
        {
          while (!get_pixel(image, bgColor, i+(height-1-j), j, THRESHOLD_BOND) && j > 0 && i+(height-1-j) < width)
            j--;
          int l = 0;
          while (get_pixel(image, bgColor, i+(height-1-j), j, THRESHOLD_BOND) && j > 0 && i+(height-1-j) < width)
            {
              l++;
              j--;
            }
          if (l < max_thick)
            n[l]++;
	  lines.push_back(l);
        }
	}*/
  double max_v = 0;
  max = 1;
  sort(lines.begin(),lines.end());
  for (int l = 1; l < max_thick; l++)
    {
      //cout << l << " " << n[l] << endl;
      if (n[l] > max_v)
        {
          max_v = n[l];
          max = l;
        }
    }
  //if (lines.size() > 1)
  // max = lines[lines.size() / 2];

  if (max > 2)
    nf = n[2] / n[3];
  else if (max == 2)
    nf = n[1] / n[2];
  else
    nf = n[2] / n[1];
  if (n[5]!=0) nf45=n[4]/n[5];
  else nf45=0;

  return (nf);
}
