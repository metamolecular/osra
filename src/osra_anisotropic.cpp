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

#define cimg_use_magick
#define cimg_plugin "greycstoration.h"

#include "CImg.h"

using namespace cimg_library;
using namespace Magick;

Image anisotropic_smoothing(const Image &image, int width, int height, const float amplitude, const float sharpness,
                            const float anisotropy, const float alpha, const float sigma)
{

  Image res(Geometry(width, height), "white");
  res.type(GrayscaleType);
  #pragma omp critical
  {
    CImg<unsigned char> source(width, height, 1, 1, 0);
    unsigned char color[1] = { 0 };
    unsigned char cc;
    ColorGray c;

    for (int i = 0; i < width; i++)
      for (int j = 0; j < height; j++)
        {
          c = image.pixelColor(i, j);
          color[0] = (unsigned char) (255 * c.shade());
          source.draw_point(i, j, color);
        }
    CImg<unsigned char> dest(source);
    const float gfact = 1.;
    //const float amplitude = 5.; // 20
    // const float sharpness = 0.3;
    //const float anisotropy = 1.;
    //const float alpha = .2; //0.6
    //const float sigma = 1.1; // 2.
    const float dl = 0.8;
    const float da = 30.;
    const float gauss_prec = 2.;
    const unsigned int interp = 0;
    const bool fast_approx = true;
    const unsigned int tile = 512;
    const unsigned int btile = 4;
    const unsigned int threads = 1; // orig - 2

    dest.greycstoration_run(amplitude, sharpness, anisotropy, alpha, sigma, gfact, dl, da, gauss_prec, interp,
                            fast_approx, tile, btile, threads);
    do
      {
        cimg::wait(1);
      }
    while (dest.greycstoration_is_running());

    for (int i = 0; i < width; i++)
      for (int j = 0; j < height; j++)
        {
          cc = dest(i, j);
          c.shade(1. * cc / 255);
          res.pixelColor(i, j, c);
        }
  }
  return (res);
}

Image anisotropic_scaling(const Image &image, int width, int height, int nw, int nh)
{

  Image res(Geometry(nw, nh), "white");
  res.type(GrayscaleType);
  #pragma omp critical
  {
    CImg<unsigned char> source(width, height, 1, 1, 0);
    unsigned char color[1] = { 0 };
    unsigned char cc;
    ColorGray c;


    for (int i = 0; i < width; i++)
      for (int j = 0; j < height; j++)
        {
          c = image.pixelColor(i, j);
          color[0] = (unsigned char) (255 * c.shade());
          source.draw_point(i, j, color);
        }

    //const float gfact = (sizeof(T) == 2) ? 1.0f / 256 : 1.0f;
    const float gfact = 1.;
    const float amplitude = 20.; // 40 20!
    const float sharpness = 0.2; // 0.2! 0.3
    const float anisotropy = 1.;
    const float alpha = .6; //0.6! 0.8
    const float sigma = 2.; //1.1 2.!
    const float dl = 0.8;
    const float da = 30.;
    const float gauss_prec = 2.;
    const unsigned int interp = 0;
    const bool fast_approx = true;
    const unsigned int tile = 512; // 512 0
    const unsigned int btile = 4;
    const unsigned int threads = 1; // 2 1

    const unsigned int init = 5;
    CImg<unsigned char> mask;

    mask.assign(source.dimx(), source.dimy(), 1, 1, 255);
    mask = !mask.resize(nw, nh, 1, 1, 4);
    source.resize(nw, nh, 1, -100, init);
    CImg<unsigned char> dest(source);

    dest.greycstoration_run(mask, amplitude, sharpness, anisotropy, alpha, sigma, gfact, dl, da, gauss_prec, interp,
                            fast_approx, tile, btile, threads);
    do
      {
        cimg::wait(1);
      }
    while (dest.greycstoration_is_running());

    for (int i = 0; i < nw; i++)
      for (int j = 0; j < nh; j++)
        {
          cc = dest(i, j);
          c.shade(1. * cc / 255);
          res.pixelColor(i, j, c);
        }
  }
  return (res);
}
