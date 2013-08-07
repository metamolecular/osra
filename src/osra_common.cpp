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

// File: osra_common.cpp
//
// Common functionality routines
//

#include <math.h> // fabs(double)
#include <float.h> // FLT_MAX
#include <fstream> // std::ofstream, std::ifstream

#include "osra_segment.h"
#include "osra_common.h"

/* return new un-initialized bitmap. NULL with errno on error */
potrace_bitmap_t *const bm_new(int w, int h)
{
  potrace_bitmap_t *bm;
  int dy = (w + BM_WORDBITS - 1) / BM_WORDBITS;

  bm = (potrace_bitmap_t *) malloc(sizeof(potrace_bitmap_t));
  if (!bm)
    {
      return NULL;
    }
  bm->w = w;
  bm->h = h;
  bm->dy = dy;
  bm->map = (potrace_word *) malloc(dy * h * BM_WORDSIZE);
  if (!bm->map)
    {
      free(bm);
      return NULL;
    }
  return bm;
}

double distance(double x1, double y1, double x2, double y2)
{
  return (sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
}

double angle4(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
  double p, l1, l2, cos;

  p = (x1 - x2) * (x3 - x4) + (y1 - y2) * (y3 - y4);
  l1 = distance(x1, y1, x2, y2);
  l2 = distance(x4, y4, x3, y3);
  cos = p / (l1 * l2);
  return (cos);
}

int get_pixel(const Image &image, const ColorGray &bg, unsigned int x, unsigned int y, double THRESHOLD)
{
  if ((x < image.columns()) && (y < image.rows()))
    {
      ColorGray c = image.pixelColor(x, y);
      if (fabs(c.shade() - bg.shade()) > THRESHOLD)
        return (1);
    }
  return (0);
}

void delete_curve(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond,
                  const potrace_path_t * const curve)
{
  for (int i = 0; i < n_atom; i++)
    {
      if (atom[i].curve == curve)
        {
          atom[i].exists = false;
        }
    }
  for (int i = 0; i < n_bond; i++)
    {
      if (bond[i].curve == curve)
        {
          bond[i].exists = false;
        }
    }
}

void delete_curve_with_children(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond,
                                const potrace_path_t * const p)
{
  delete_curve(atom, bond, n_atom, n_bond, p);
  potrace_path_t *child = p->childlist;
  while (child != NULL)
    {
      delete_curve(atom, bond, n_atom, n_bond, child);
      child = child->sibling;
    }
}

double angle_between_bonds(const vector<bond_t> &bond, int i, int j, const vector<atom_t> &atom)
{
  return (angle4(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[j].a].x,
                 atom[bond[j].a].y, atom[bond[j].b].x, atom[bond[j].b].y));
}

double bond_length(const vector<bond_t> &bond, int i, const vector<atom_t> &atom)
{
  return (distance(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x, atom[bond[i].b].y));
}

double distance_from_bond_y(double x0, double y0, double x1, double y1, double x, double y)
{
  double d1 = distance(x0, y0, x1, y1);
  double cos = (x1 - x0) / d1;
  double sin = -(y1 - y0) / d1;
  double h = -(x - x0) * sin - (y - y0) * cos;
  return (h);
}

double distance_between_bonds(const vector<bond_t> &bond, int i, int j, const vector<atom_t> &atom)
{
  /*
  double y1 = distance_from_bond_y(atom[bond[j].a].x, atom[bond[j].a].y, atom[bond[j].b].x, atom[bond[j].b].y,
  		atom[bond[i].a].x, atom[bond[i].a].y);
  double y2 = distance_from_bond_y(atom[bond[j].a].x, atom[bond[j].a].y, atom[bond[j].b].x, atom[bond[j].b].y,
  		atom[bond[i].b].x, atom[bond[i].b].y);
  if (fabs(y1 - y2) >= 4)
  	return (FLT_MAX);
  double r1 = max(fabs(y1), fabs(y2));
  */
  double y3 = distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x, atom[bond[i].b].y,
                                   atom[bond[j].a].x, atom[bond[j].a].y);
  double y4 = distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x, atom[bond[i].b].y,
                                   atom[bond[j].b].x, atom[bond[j].b].y);
  if (fabs(y3 - y4) >= 4)
    return (FLT_MAX);
  double r2 = max(fabs(y3), fabs(y4));
  return (r2);
}

double distance_from_bond_x_a(double x0, double y0, double x1, double y1, double x, double y)
{
  double d1 = distance(x0, y0, x1, y1);
  double cos = (x1 - x0) / d1;
  double sin = -(y1 - y0) / d1;
  double l = (x - x0) * cos - (y - y0) * sin;
  return (l);
}

double distance_from_bond_x_b(double x0, double y0, double x1, double y1, double x, double y)
{
  double d1 = distance(x0, y0, x1, y1);
  double cos = (x1 - x0) / d1;
  double sin = -(y1 - y0) / d1;
  double l = (x - x0) * cos - (y - y0) * sin;
  return (l - d1);
}

double percentile75(const vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom)
{
  vector<double> a;
  int n = 0;

  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      {
        a.push_back(bond_length(bond, i, atom));
        n++;
      }
  if (n > 1)
    {
      std::sort(a.begin(), a.end());
      int pos = 3 * (n - 1) / 4;
      return (a[pos]);
    }
  else
    return (10.0);
}


bool terminal_bond(int a, int b, const vector<bond_t> &bond, int n_bond)
{
  bool terminal = true;

  for (int l = 0; l < n_bond; l++)
    if (l != b && bond[l].exists && (bond[l].a == a || bond[l].b == a))
      terminal = false;

  return (terminal);
}


void debug_image(Image image, const vector<atom_t> &atom, int n_atom, const vector<bond_t> &bond, int n_bond,
               const string &fname)
{
  image.modifyImage();
  image.type(TrueColorType);
  image.strokeWidth(1);

  int max_x = image.columns();
  int max_y = image.rows();

  for (int i = 0; i < n_bond; i++)
    {
      if ((bond[i].exists) && (atom[bond[i].a].exists) && (atom[bond[i].b].exists))
        {
          if (bond[i].type == 1)
            {
              image.strokeColor("green");
            }
          else if (bond[i].type == 2)
            {
              image.strokeColor("yellow");
            }
          else if (bond[i].type >= 3)
            {
              image.strokeColor("red");
            }
          if (bond[i].hash)
            {
              image.strokeColor("blue");
            }
          else if (bond[i].wedge)
            {
              image.strokeColor("purple");
            }
          image.draw(DrawableLine(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x, atom[bond[i].b].y));
        }
    }

  for (int i = 0; i < n_atom; i++)
    {
      if (atom[i].exists)
        {
          if ((int(atom[i].x) < max_x) && (int(atom[i].y < max_y)))
            image.pixelColor(int(atom[i].x), int(atom[i].y), "blue");
        }
    }

  image.write(fname);
}

void draw_square(Image &image, int x1, int y1, int x2, int y2, const string &color)
{
  image.strokeWidth(1);
  image.strokeColor(color);
  image.draw(DrawableLine(x1, y1, x2, y1));
  image.draw(DrawableLine(x1, y2, x2, y2));
  image.draw(DrawableLine(x1, y1, x1, y2));
  image.draw(DrawableLine(x2, y1, x2, y2));
}

void draw_box(Image &image, vector<box_t> &boxes, int n_boxes, const string &fname)
{
  image.modifyImage();
  image.type(TrueColorType);

  for (int i = 0; i < n_boxes; i++)
    {
      draw_square(image, boxes[i].x1, boxes[i].y1, boxes[i].x2, boxes[i].y2, "green");
    }
  image.write(fname);
}


int count_pages(const string &input)
{
  list<Image> imageList;
  readImages(&imageList, input);
  return (imageList.size());
}

int count_atoms(const vector<atom_t> &atom, int n_atom)
{
  int r = 0;
  for (int i = 0; i < n_atom; i++)
    if (atom[i].exists)
      r++;
  return (r);
}

int count_bonds(const vector<bond_t> &bond, int n_bond, int &bond_max_type)
{
  int r = 0;
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      {
        r++;
        if (bond[i].type>bond_max_type) bond_max_type = bond[i].type;
      }
  return (r);
}

bool detect_curve(vector<bond_t> &bond, int n_bond, const potrace_path_t * const curve)
{
  bool res = false;
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && bond[i].curve == curve && bond[i].type == 1 && !bond[i].wedge && !bond[i].hash)
      res = true;
  return (res);
}


// Igor Filippov - 2009.
// The following two functions are adapted from ConfigFile
///
// Class for reading named values from configuration files
// Richard J. Wagner  v2.1  24 May 2004  wagnerr@umich.edu

// Copyright (c) 2004 Richard J. Wagner
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

void trim(string &s)
{
  // Remove leading and trailing whitespace
  static const char whitespace[] = " \n\t\v\r\f";
  s.erase(0, s.find_first_not_of(whitespace));
  s.erase(s.find_last_not_of(whitespace) + 1U);
}

bool load_config_map(const string &file, map<string, string> &out)
{
  typedef string::size_type pos;
  const string& delim = " "; // separator
  const pos skip = delim.length(); // length of separator

  std::ifstream is(file.c_str());
  if (!is)
    return false;

  while (is)
    {
      // Read an entire line at a time
      string line;
      std::getline(is, line);

      // Ignore comments
      //line = line.substr(0, line.find(comm));
      if (line.length() == 0 || line.at(0) == '#')
        continue;

      // replace tabs with spaces
      pos t;
      while ((t = line.find('\t')) != string::npos)
        line[t] = ' ';

      // Parse the line if it contains a delimiter
      pos delimPos = line.find(delim);
      if (delimPos < string::npos)
        {
          // Extract the key
          string key = line.substr(0, delimPos);
          line.replace(0, delimPos + skip, "");

          // Store key and value
          trim(key);
          trim(line);
          out[key] = line; // overwrites if key is repeated
        }
    }

  is.close();

  return true;
}


bool comp_boxes(const box_t &aa, const box_t &bb)
{
  if (aa.y2 < bb.y1)
    return (true);
  if (aa.y1 > bb.y2)
    return (false);
  if (aa.x1 > bb.x1)
    return (false);
  if (aa.x1 < bb.x1)
    return (true);
  return (false);
}

