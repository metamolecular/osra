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

// File: osra_grayscale.cpp
//
// Defines grayscale conversion functions
//

#include <vector>
#include <iostream> // std::ostream, std::cout
#include <math.h> // fabs(double)

#include "osra.h"
#include "osra_grayscale.h"

const Color getBgColor(const Image &image)
{
  ColorGray c, r;
  r = image.pixelColor(1, 1);
  for (int i = 0; i < BG_PICK_POINTS; i++)
    {
      double a = (double)rand()/RAND_MAX;
      double b = (double)rand()/RAND_MAX;
      int x = int(image.columns() *a);
      int y = int(image.rows() * b);
      c = image.pixelColor(x, y);
      if (c.shade() > r.shade())
        r = c;
    }

  return (r);
}

void otsu_find_peaks(const vector<int> &h, int num_bins, int &peak1, int &peak2, int &max1, int &max2)
{
// Otsu Algorithm, from http://habrahabr.ru/blogs/algorithm/112079/
  unsigned int m = 0;
  unsigned int n = 0;
  double maxSigma = -1;
  unsigned int min_t = num_bins;
  unsigned int alpha1 = 0;
  unsigned int beta1 = 0;
  for (unsigned int i = 0; i<num_bins; i++)
    {
      m += i*h[i];
      n += h[i];
    }

  for (unsigned int t = 0; t < num_bins; t++)
    {
      alpha1 += t * h[t];
      beta1 += h[t];

      double w1 = (double)beta1 / n;
      double a = (double)alpha1 / beta1 - (double)(m - alpha1) / (n - beta1);
      double sigma = w1 * (1 - w1) * a * a;

      if (sigma > maxSigma)
        {
          maxSigma = sigma;
          min_t = t;
        }
    }
  max1 = 0;
  peak1 = 0;
  for (unsigned int i = 0; i<=min_t; i++)
    if (h[i] > max1)
      {
        max1 = h[i];
        peak1 = i;
      }
  max2 = 0;
  peak2 = 0;
  for (unsigned int i = min_t+1; i<num_bins; i++)
    if (h[i] > max2)
      {
        max2 = h[i];
        peak2 = i;
      }
}

Image adaptive_otsu(const Image &image, int window)
{
  int num_bins=20;
  Image result(Geometry(image.columns(),image.rows()),"white");
  vector<int> h(num_bins,0);
  vector<int> h0(num_bins,0);
  ColorGray g;
  int peak1, peak2, max1, max2;

  for (int i1 = 0; i1 < min((int)image.columns(), window/2); i1++)
    for (int j1 = 0; j1 < min((int)image.rows(), window/2); j1++)
      {
        g = image.pixelColor(i1, j1);
        h0[int((num_bins-1)*g.shade())]++;
      }

  for (int j = 0; j < image.rows(); j++)
    {
      for (int k = 0; k<num_bins; k++) h[k]=h0[k];
      for (int i = 0; i < image.columns(); i++)
        {
          otsu_find_peaks(h,num_bins,peak1,peak2, max1, max2);
          g = image.pixelColor(i, j);
          double median = 0.5*(peak1+peak2)/num_bins;
          if (g.shade() > median)
            result.pixelColor(i,j,"white");
          else
            result.pixelColor(i,j,"black");
          if ((i-window/2) >=0)
            for (int j1 = max(0,j-window/2); j1 < min((int)image.rows(), j + window/2); j1++)
              {
                g = image.pixelColor(i-window/2, j1);
                h[int((num_bins-1)*g.shade())]--;
              }
          if ((i+window/2) < image.columns())
            for (int j1 = max(0,j-window/2); j1 < min((int)image.rows(), j + window/2); j1++)
              {
                g = image.pixelColor(i+window/2, j1);
                h[int((num_bins-1)*g.shade())]++;
              }
        }
      if ((j-window/2) >=0)
        for (int i1 = 0; i1 < min((int)image.columns(),window/2); i1++)
          {
            g = image.pixelColor(i1,j-window/2);
            h0[int((num_bins-1)*g.shade())]--;
          }
      if ((j+window/2) < image.rows())
        for (int i1 = 0; i1 < min((int)image.columns(),window/2); i1++)
          {
            g = image.pixelColor(i1,j+window/2);
            h0[int((num_bins-1)*g.shade())]++;
          }
    }
  return(result);
}

bool convert_to_gray(Image &image, bool invert, bool adaptive, bool verbose)
{
  int num_bins=50;
  int num_bins_rgb = 20;
  vector<int> h(num_bins,0);
  vector < vector < vector <int> > > bg_search(num_bins_rgb, vector < vector <int> > (num_bins_rgb, vector<int>(num_bins_rgb, 0)));
  ColorRGB c,b;
  Color t;
  ColorGray g;
  double a;

  image.type(TrueColorMatteType);
  for (int i = 0; i < BG_PICK_POINTS; i++)
    {
      double a = (double) rand() / RAND_MAX;
      double b = (double) rand() / RAND_MAX;
      int x = int(image.columns() * a);
      int y = int(image.rows() * b);
      c = image.pixelColor(x, y);
      bg_search[int((num_bins_rgb-1)*c.red())][int((num_bins_rgb-1)*c.green())][int((num_bins_rgb-1)*c.blue())]++;
    }
  int bg_peak = 0;
  double bg_pos_red = 0, bg_pos_green = 0, bg_pos_blue = 0;
  for (int i=0; i<num_bins_rgb; i++)
    for (int j=0; j<num_bins_rgb; j++)
      for (int k=0; k<num_bins_rgb; k++)
        if (bg_search[i][j][k] > bg_peak)
          {
            bg_peak = bg_search[i][j][k];
            bg_pos_red = (double)i/(num_bins_rgb-1);
            bg_pos_green = (double)j/(num_bins_rgb-1);
            bg_pos_blue = (double)k/(num_bins_rgb-1);
          }

  bool color_background = false;
  if (verbose)
    {
      cout<<"Background rgb: "<<bg_pos_red<<" "<<bg_pos_green<<" "<<bg_pos_blue<<endl;
    }

  if (fabs(bg_pos_red-bg_pos_green) > 0.05 || fabs(bg_pos_red-bg_pos_blue)>0.05 || fabs(bg_pos_green-bg_pos_blue)>0.05) color_background = true;

  bool matte = image.matte();
  if (color_background)
    {
      image.contrast(2);
      image.type(GrayscaleType);
    }

  for (unsigned int i = 0; i < image.columns(); i++)
    for (unsigned int j = 0; j < image.rows(); j++)
      {
        t = image.pixelColor(i, j);
        b = t;
        g = t;
        if (matte && t.alpha() == 1 && g.shade() < 0.5)
          {
            g.shade(1);
            image.pixelColor(i, j, g);
          }
        else if (!color_background &&
                 (fabs(b.red()-b.green()) > 0.1 || fabs(b.red()-b.blue()) > 0.1  || fabs(b.blue()-b.green()) > 0.1))
          {
            if (fabs(b.red()-bg_pos_red) >= fabs(b.green()-bg_pos_green) && fabs(b.red()-bg_pos_red) >= fabs(b.blue()-bg_pos_blue))
              a = b.red();
            else if (fabs(b.red()-bg_pos_red) < fabs(b.green()-bg_pos_green) && fabs(b.green()-bg_pos_green) >= fabs(b.blue()-bg_pos_blue))
              a = b.green();
            else
              a = b.blue();
            c.red(a);
            c.green(a);
            c.blue(a);
            image.pixelColor(i, j, c);
          }
        g = image.pixelColor(i, j);
        h[int((num_bins-1)*g.shade())]++;
      }

  int peak1, peak2, max1, max2;
  otsu_find_peaks(h,num_bins,peak1,peak2, max1, max2);

  double distance_between_peaks = (double)(peak2-peak1)/(num_bins-1);
  //  if (distance_between_peaks < THRESHOLD_GLOBAL) adaptive = true;
  if (distance_between_peaks < 0.5) adaptive = true;

  if (max1 > max2 || invert)
    invert = true;

  if (verbose)
    {
      cout << "Distance between light and dark: " << distance_between_peaks << endl;
      cout<<"Max at peak 1: "<<max1<<"  Max at peak 2: "<<max2<<endl;
      cout<<"Color background? "<<color_background<<endl;
      cout<<"Adaptive? "<<adaptive<<endl;
      cout<<"Invert? "<<invert<<endl;
    }

  //const double kernel[]={0.0, -1.0, 0.0,-1.0, 5.0, -1.0, 0.0, -1.0, 0.0};
  //image.convolve(3,kernel);

  if (!color_background)
    {
      image.contrast(2);
      image.type(GrayscaleType);
    }

  int window = min(image.columns(),image.rows()) / 41;
  if (window < 15) window = 15;

  if (adaptive)
    {
      image.despeckle();
      if (invert)
        {
          image.adaptiveThreshold(window,window,7);
        }
      else
        {
          image.negate();
          image.adaptiveThreshold(window,window,7);
          image.negate();
        }
    }
  else if (color_background)
    {
      image.despeckle();
      image = adaptive_otsu(image,window);
    }

  if (invert)
    image.negate();

  //  image.write("tmp.png");

  return(adaptive);
}
