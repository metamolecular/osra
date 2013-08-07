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

// File: osra_labels.cpp
//
// Defines functions dealing with atomic labels
//
#include <math.h> // fabs(double)
#include <float.h> // FLT_MAX
#include <limits.h> // INT_MAX
#include <iostream> // std::ostream, std::cout

#include "osra_common.h"
#include "osra_ocr.h"
#include "osra_labels.h"


bool alone(const vector<bond_t> &bond, int i, double avg)
{
  bool alone = false;
  const potrace_path_t * const p = bond[i].curve;

  if ((p->sign == int('+')) && (p->area < 2 * avg))
    alone = true;

  return (alone);
}

void delete_bonds_in_char(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double left, double top,
                          double right, double bottom)
{
  for (int j = 0; j < n_bond; j++)
    if (bond[j].exists && atom[bond[j].a].x >= left && atom[bond[j].a].x <= right && atom[bond[j].a].y >= top
        && atom[bond[j].a].y <= bottom && atom[bond[j].b].x >= left && atom[bond[j].b].x <= right
        && atom[bond[j].b].y >= top && atom[bond[j].b].y <= bottom)
      bond[j].exists = false;
}

bool chlorine(const vector<bond_t> &bond, const vector<atom_t> &atom, int i, vector<letters_t> &letters, int n_letters,
              int max_font_height, int min_font_height)
{
  bool res = false;
  double x = (atom[bond[i].a].x + atom[bond[i].b].x) / 2;
  double y = (atom[bond[i].a].y + atom[bond[i].b].y) / 2;
  double r = bond_length(bond, i, atom) / 2;
  if ((bond_length(bond, i, atom) < max_font_height) && (bond_length(bond, i, atom) > min_font_height) && (fabs(
        atom[bond[i].a].x - atom[bond[i].b].x) < fabs(atom[bond[i].a].y - atom[bond[i].b].y)))
    {
      for (int j = 0; j < n_letters; j++)
        {
          if ((distance(x, y, letters[j].x, letters[j].y) < r + letters[j].r) && (fabs(y - letters[j].y) < min(r,letters[j].r)))
            {
              res = true;
            }
        }
    }

  return (res);
}

bool iodine(const vector<bond_t> &bond, const vector<atom_t> &atom, int i, vector<letters_t> &letters, int n_letters, int max_font_height, int min_font_height)
{
  if (bond_length(bond, i, atom) < max_font_height && 
      bond_length(bond, i, atom) > min_font_height &&
      fabs(atom[bond[i].a].x - atom[bond[i].b].x) < V_DISPLACEMENT && 
      fabs(atom[bond[i].a].y - atom[bond[i].b].y) > min_font_height &&
      bond[i].type == 1)
    {
      bool found = false;
      double xa = atom[bond[i].a].x;
      double ya = atom[bond[i].a].y;
      double xb = atom[bond[i].b].x;
      double yb = atom[bond[i].b].y;
      double bl = bond_length(bond, i, atom);

      for (int j = 0; j < n_letters; j++)
	{
	  double d1 = fabs(distance_from_bond_x_a(xa, ya, xb, yb, letters[j].x, letters[j].y));
	  double d2 = fabs(distance_from_bond_x_b(xa, ya, xb, yb, letters[j].x, letters[j].y));
	  double h = fabs(distance_from_bond_y(xa, ya, xb, yb, letters[j].x, letters[j].y));
	  double nb = min(d1,d2) - letters[j].r;

	  if (nb <= bl/2 && h <=  letters[j].r / 2)
	    found = true;
	}
      if (!found)
	return true;
    }

  return false;
}

int remove_small_bonds(vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom, vector<letters_t> &letters,
                       int n_letters, int max_font_height, int min_font_height, double avg)
{
  for (int i = 0; i < n_bond; i++)
    if ((bond[i].exists) && (bond[i].type == 1))
      {
        bool al = alone(bond, i, avg);
        if (bond_length(bond, i, atom) < V_DISPLACEMENT)
          {
            bond[i].exists = false;
          }
        else if (al) 
          {
	    bool cl = chlorine(bond, atom, i, letters, n_letters, max_font_height, min_font_height);
	    bool io = iodine(bond, atom, i, letters, n_letters, max_font_height, min_font_height);
	    //if (io)
	    //cout<<atom[bond[i].a].x<<" "<<atom[bond[i].a].y<<" "<<atom[bond[i].b].x<<" "<<atom[bond[i].b].y<<endl;
	    if (cl || io)
	      {
		letters_t lt;
		letters.push_back(lt);
		if (cl)
		  letters[n_letters].a = 'l';  // small L for chlorine
		else
		  letters[n_letters].a = 'I';   // capital i for iodine
		letters[n_letters].x = (atom[bond[i].a].x + atom[bond[i].b].x) / 2;
		letters[n_letters].y = (atom[bond[i].a].y + atom[bond[i].b].y) / 2;
		letters[n_letters].r = bond_length(bond, i, atom) / 2;
		letters[n_letters].min_x = min(atom[bond[i].a].min_x,atom[bond[i].b].min_x);
		letters[n_letters].max_x = max(atom[bond[i].b].max_x,atom[bond[i].b].max_x);
		letters[n_letters].min_y = min(atom[bond[i].a].min_y,atom[bond[i].b].min_y);
		letters[n_letters].max_y = max(atom[bond[i].a].max_y,atom[bond[i].b].max_y);
		letters[n_letters].free = true;
		n_letters++;
		if (n_letters >= MAX_ATOMS)
		  n_letters--;
		bond[i].exists = false;
	      }
          }
      }
  return (n_letters);
}

bool comp_lbonds(const lbond_t &left, const lbond_t &right)
{
  if (left.x < right.x)
    return (true);
  if (left.x > right.x)
    return (false);
  return (false);
}

bool comp_letters(const letters_t &left, const letters_t &right)
{
  if (left.x < right.x)
    return (true);
  if (left.x > right.x)
    return (false);
  return (false);
}

int assemble_labels(vector<letters_t> &letters, int n_letters, vector<label_t> &label)
{
  vector<lbond_t> lbond;

  std::sort(letters.begin(), letters.end(), comp_letters);

  for (int i = 0; i < letters.size(); i++)
    {
      //            cout<<letters[i].a<<" "<<letters[i].min_x<<" "<<letters[i].min_y<<" "<<letters[i].max_x<<" "<<letters[i].max_y<<endl;
      for (int j = i + 1; j < letters.size(); j++)
	{
	  bool digits = fabs(letters[i].y - letters[j].y) < (letters[i].r + letters[j].r) 
							    && ((letters[i].y < letters[j].y  && isdigit(letters[j].a)) || ((letters[j].y < letters[i].y) && isdigit(letters[i].a)));
	  bool alphanum = distance(letters[i].x, letters[i].y, letters[j].x, letters[j].y) < 2 * max(letters[i].r, letters[j].r)
											     && (fabs(letters[i].y - letters[j].y) < min(letters[i].r, letters[j].r) 
												 || fabs(letters[i].x - letters[j].x) < min(letters[i].r, letters[j].r)
												 || digits);
	  bool signs = distance(letters[i].x, letters[i].y, letters[j].x, letters[j].y) < 1.5 * (letters[i].r + letters[j].r) 
											  && ((letters[i].y < letters[j].y && (letters[i].a == '-' || letters[i].a == '+')) || 
											      (letters[j].y < letters[i].y && (letters[j].a == '-' || letters[j].a == '+')));
	  if (alphanum || signs)
	    {
	      lbond_t lb;
	      lb.a = i;
	      lb.b = j;
	      lb.x = letters[i].x;
	      lb.exists = true;
	      lbond.push_back(lb);
	      letters[i].free = false;
	      letters[j].free = false;
	      break;
	    }
	}
    }

  std::sort(lbond.begin(), lbond.end(), comp_lbonds);

  for (int i = 0; i < lbond.size(); i++)
    if (lbond[i].exists)
      {
        bool found_left = false;
        label_t lb;

        lb.x1 = FLT_MAX;
        lb.y1 = FLT_MAX;
        lb.r1 = 0;
        lb.x2 = FLT_MAX;
        lb.y2 = FLT_MAX;
        lb.r2 = 0;
        lb.a = letters[lbond[i].a].a;
        lb.a += letters[lbond[i].b].a;
        lb.n.push_back(lbond[i].a);
        lb.n.push_back(lbond[i].b);
	lb.min_x =  min(letters[lbond[i].a].min_x, letters[lbond[i].b].min_x);
	lb.min_y =  min(letters[lbond[i].a].min_y,letters[lbond[i].b].min_y);
	lb.max_x =  max(letters[lbond[i].a].max_x, letters[lbond[i].b].max_x);
	lb.max_y =  max(letters[lbond[i].a].max_y,letters[lbond[i].b].max_y);

        if (!isdigit(letters[lbond[i].a].a) && letters[lbond[i].a].a != '-' && letters[lbond[i].a].a != '+'
            && !found_left)
          {
            lb.x1 = letters[lbond[i].a].x;
            lb.y1 = letters[lbond[i].a].y;
            lb.r1 = letters[lbond[i].a].r;
            found_left = true;
          }
        if (!isdigit(letters[lbond[i].b].a) && letters[lbond[i].b].a != '-' && letters[lbond[i].b].a != '+'
            && !found_left)
          {
            lb.x1 = letters[lbond[i].b].x;
            lb.y1 = letters[lbond[i].b].y;
            lb.r1 = letters[lbond[i].b].r;
            found_left = true;
          }
        if (!isdigit(letters[lbond[i].a].a) && letters[lbond[i].a].a != '-' && letters[lbond[i].a].a != '+')
          {
            lb.x2 = letters[lbond[i].a].x;
            lb.y2 = letters[lbond[i].a].y;
            lb.r2 = letters[lbond[i].a].r;
          }
        if (!isdigit(letters[lbond[i].b].a) && letters[lbond[i].b].a != '-' && letters[lbond[i].b].a != '+')
          {
            lb.x2 = letters[lbond[i].b].x;
            lb.y2 = letters[lbond[i].b].y;
            lb.r2 = letters[lbond[i].b].r;
          }
        lbond[i].exists = false;
        int last = lbond[i].b;
        for (int j = i + 1; j < lbond.size(); j++)
          if ((lbond[j].exists) && (lbond[j].a == last))
            {
              lb.a += letters[lbond[j].b].a;
	      lb.min_x =  min(lb.min_x, letters[lbond[j].b].min_x);
	      lb.min_y =  min(lb.min_y,letters[lbond[j].b].min_y);
	      lb.max_x =  max(lb.max_x, letters[lbond[j].b].max_x);
	      lb.max_y =  max(lb.max_y,letters[lbond[j].b].max_y);
              lb.n.push_back(lbond[j].b);
              if (!isdigit(letters[lbond[j].a].a) && letters[lbond[j].a].a != '-' && letters[lbond[j].a].a != '+'
                  && !found_left)
                {
                  lb.x1 = letters[lbond[j].a].x;
                  lb.y1 = letters[lbond[j].a].y;
                  lb.r1 = letters[lbond[j].a].r;
                  found_left = true;
                }
              if (!isdigit(letters[lbond[j].b].a) && letters[lbond[j].b].a != '-' && letters[lbond[j].b].a != '+'
                  && !found_left)
                {
                  lb.x1 = letters[lbond[j].b].x;
                  lb.y1 = letters[lbond[j].b].y;
                  lb.r1 = letters[lbond[j].b].r;
                  found_left = true;
                }
              if (!isdigit(letters[lbond[j].a].a) && letters[lbond[j].a].a != '-' && letters[lbond[j].a].a != '+')
                {
                  lb.x2 = letters[lbond[j].a].x;
                  lb.y2 = letters[lbond[j].a].y;
                  lb.r2 = letters[lbond[j].a].r;
                }
              if (!isdigit(letters[lbond[j].b].a) && letters[lbond[j].b].a != '-' && letters[lbond[j].b].a != '+')
                {
                  lb.x2 = letters[lbond[j].b].x;
                  lb.y2 = letters[lbond[j].b].y;
                  lb.r2 = letters[lbond[j].b].r;
                }
              last = lbond[j].b;
              lbond[j].exists = false;
            }
        label.push_back(lb);
      }

   int old_n_label = label.size();
  for (int i = 0; i < old_n_label; i++)
    {
      double cy = 0;
      int n = 0;

      for (unsigned int j = 0; j < label[i].n.size(); j++)
        if (isalpha(letters[label[i].n[j]].a))
          {
            cy += letters[label[i].n[j]].y;
            n++;
          }
      cy /= n;
      n = 0;
      for (unsigned int j = 0; j < label[i].n.size(); j++)
        if (isalpha(letters[label[i].n[j]].a) && letters[label[i].n[j]].y - cy > letters[label[i].n[j]].r / 2)
          n++;

      if (n > 1)
        {
	  vector<int> old_label_n(label[i].n);
          label[i].a = "";
          label[i].x1 = FLT_MAX;
          label[i].x2 = 0;
	  label[i].n.clear();
          label_t lb;

          lb.a = "";
          lb.x1 = FLT_MAX;
          lb.x2 = 0;
	  lb.min_x = INT_MAX;
	  lb.min_y = INT_MAX;
	  lb.max_x = 0;
	  lb.max_y = 0;
          for (unsigned int j = 0; j < old_label_n.size(); j++)
            {
              if (letters[old_label_n[j]].y > cy)
                {
                  label[i].a += letters[old_label_n[j]].a;
		  label[i].n.push_back(old_label_n[j]);
		  label[i].min_x =  min(label[i].min_x, letters[old_label_n[j]].min_x);
		  label[i].min_y =  min(label[i].min_y,letters[old_label_n[j]].min_y);
		  label[i].max_x =  max(label[i].max_x, letters[old_label_n[j]].max_x);
		  label[i].max_y =  max(label[i].max_y,letters[old_label_n[j]].max_y);
                  if (isalpha(letters[old_label_n[j]].a))
                    {
                      if (letters[old_label_n[j]].x < label[i].x1)
                        {
                          label[i].x1 = letters[old_label_n[j]].x;
                          label[i].y1 = letters[old_label_n[j]].y;
                          label[i].r1 = letters[old_label_n[j]].r;
                        }
                      if (letters[old_label_n[j]].x > label[i].x2)
                        {
                          label[i].x2 = letters[old_label_n[j]].x;
                          label[i].y2 = letters[old_label_n[j]].y;
                          label[i].r2 = letters[old_label_n[j]].r;
                        }
                    }
                }
              else
                {
                  lb.a += letters[old_label_n[j]].a;
		  lb.n.push_back(old_label_n[j]);
		  lb.min_x =  min(lb.min_x, letters[old_label_n[j]].min_x);
		  lb.min_y =  min(lb.min_y,letters[old_label_n[j]].min_y);
		  lb.max_x =  max(lb.max_x, letters[old_label_n[j]].max_x);
		  lb.max_y =  max(lb.max_y,letters[old_label_n[j]].max_y);
                  if (isalpha(letters[old_label_n[j]].a))
                    {
                      if (letters[old_label_n[j]].x < lb.x1)
                        {
                          lb.x1 = letters[old_label_n[j]].x;
                          lb.y1 = letters[old_label_n[j]].y;
                          lb.r1 = letters[old_label_n[j]].r;
                        }
                      if (letters[old_label_n[j]].x > lb.x2)
                        {
                          lb.x2 = letters[old_label_n[j]].x;
                          lb.y2 = letters[old_label_n[j]].y;
                          lb.r2 = letters[old_label_n[j]].r;
                        }
                    }
                }
            }
          label.push_back(lb);
        }
    }
  
  for (int i = 0; i < label.size(); i++)
    {
      //cout<<label[i].a<<" "<<label[i].min_x<<" "<<label[i].min_y<<" "<<label[i].max_x<<" "<<label[i].max_y<<endl;
      bool cont = true;
      string charges = "";
      while (cont)
        {
          cont = false;
          string::size_type pos = label[i].a.find_first_of('-');
          if (pos != string::npos)
            {
              label[i].a.erase(pos, 1);
              charges += "-";
              cont = true;
            }
          pos = label[i].a.find_first_of('+');
          if (pos != string::npos)
            {
              label[i].a.erase(pos, 1);
              charges += "+";
              cont = true;
            }
        }
      label[i].a += charges;
      }
 
  return (label.size());
}


int find_numbers(const potrace_path_t * p, const Image &orig, vector<letters_t> &letters, vector<atom_t> &atom, vector<bond_t> &bond, 
		 int n_atom, int n_bond, int height, int width, ColorGray &bgColor, double THRESHOLD, int n_letters)
{
  int max_font_width, max_font_height;
  vector<int> widths, heights;
  for (int i=0; i<n_letters; i++)
    if (letters[i].a >= '2' && letters[i].a <= '9')
      {
	widths.push_back(letters[i].max_x - letters[i].min_x+1);
	heights.push_back(letters[i].max_y - letters[i].min_y+1);
      }
  if (widths.size() < 3) return n_letters;
  max_font_width = min(MAX_FONT_WIDTH,widths[widths.size() / 2]);
  max_font_height = min(MAX_FONT_HEIGHT,heights[heights.size() / 2]);


  int n, *tag;
  potrace_dpoint_t (*c)[3];

  while (p != NULL)
    {
      if ((p->sign == int('+')))
        {
          n = p->curve.n;
          tag = p->curve.tag;
          c = p->curve.c;
          int top = height;
          int x1 = 0;
          int left = width;
          int y1 = 0;
          int bottom = 0;
          int x2 = 0;
          int right = 0;
          int y2 = 0;
          int cx,cy;
          for (int i = 0; i < n; i++)
            {
              switch (tag[i])
                {
                case POTRACE_CORNER:
                  cx = c[i][1].x;
                  cy = c[i][1].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;

                  if (cx < left)
                    {
                      left = cx;
                      y1 = cy;
                    }
                  if (cx > right)
                    {
                      right = cx;
                      y2 = cy;
                    }
                  if (cy < top)
                    {
                      top = cy;
                      x1 = cx;
                    }
                  if (cy > bottom)
                    {
                      bottom = cy;
                      x2 = cx;
                    }
                  break;
                case POTRACE_CURVETO:
                  cx = c[i][0].x;
                  cy = c[i][0].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  if (cx < left)
                    {
                      left = cx;
                      y1 = cy;
                    }
                  if (cx > right)
                    {
                      right = cx;
                      y2 = cy;
                    }
                  if (cy < top)
                    {
                      top = cy;
                      x1 = cx;
                    }
                  if (cy > bottom)
                    {
                      bottom = cy;
                      x2 = cx;
                    }
                  cx = c[i][1].x;
                  cy = c[i][1].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  if (cx < left)
                    {
                      left = cx;
                      y1 = cy;
                    }
                  if (cx > right)
                    {
                      right = cx;
                      y2 = cy;
                    }
                  if (cy < top)
                    {
                      top = cy;
                      x1 = cx;
                    }
                  if (cy > bottom)
                    {
                      bottom = cy;
                      x2 = cx;
                    }
                  break;
                }
              cx = c[i][2].x;
              cy = c[i][2].y;
              if (cx<0) cx=0;
              if (cx>width) cx=width;
              if (cy<0) cy=0;
              if (cy>height) cy=height;
              if (cx < left)
                {
                  left = cx;
                  y1 = cy;
                }
              if (cx > right)
                {
                  right = cx;
                  y2 = cy;
                }
              if (cy < top)
                {
                  top = cy;
                  x1 = cx;
                }
              if (cy > bottom)
                {
                  bottom = cy;
                  x2 = cx;
                }
            }

          if (((bottom - top) <= max_font_height) && ((right - left) <= max_font_width) && (right - left
              > V_DISPLACEMENT) && (bottom - top > MIN_FONT_HEIGHT))
            {
              int s = 1;
              while ((top > 0) && (s > 0))
                {
                  s = 0;
                  s = get_pixel(orig, bgColor, x1, top, THRESHOLD);
                  if (s > 0)
                    top--;
                }
              s = 1;
              while ((bottom < height) && (s > 0))
                {
                  s = 0;
                  s = get_pixel(orig, bgColor, x2, bottom, THRESHOLD);
                  if (s > 0)
                    bottom++;
                }
              s = 1;
              while ((left > 0) && (s > 0))
                {
                  s = 0;
                  s = get_pixel(orig, bgColor, left, y1, THRESHOLD);
                  if (s > 0)
                    left--;
                }
              s = 1;
              while ((right < width) && (s > 0))
                {
                  s = 0;
                  s = get_pixel(orig, bgColor, right, y2, THRESHOLD);
                  if (s > 0)
                    right++;
                }
            }

          if (((bottom - top) <= max_font_height) && ((right - left) <= max_font_width) && (right - left
              > V_DISPLACEMENT) && (bottom - top > MIN_FONT_HEIGHT))
            {
	      bool found = false;
	      for (int i=0; i<n_letters; i++)
		{
		  int x =  (left + right) / 2;
		  int y = (top + bottom) / 2;
		  if (distance(x,y,letters[i].x,letters[i].y) < V_DISPLACEMENT)
		    {
		      found = true;
		      break;
		    }
		}
	      if (!found)
		{
		  char label = 0;
		  label = get_atom_label(orig, bgColor, left, top, right, bottom, THRESHOLD, (right + left) / 2, top, false, false ,true);
		  
		  if (label == '1')
		    {
		      letters_t lt;
		      letters.push_back(lt);
		      letters[n_letters].a = label;
		      letters[n_letters].x = (left + right) / 2;
		      letters[n_letters].y = (top + bottom) / 2;
		      letters[n_letters].r = distance(left, top, right, bottom) / 2;
		      letters[n_letters].min_x = left;
		      letters[n_letters].max_x = right;
		      letters[n_letters].min_y = top;
		      letters[n_letters].max_y = bottom;
		      letters[n_letters].free = true;
		      n_letters++;
		      if (n_letters >= MAX_ATOMS)
			n_letters--;
		      delete_bonds_in_char(bond, n_bond, atom, left, top, right, bottom);
		      delete_curve_with_children(atom, bond, n_atom, n_bond, p);
		    }
		}
            }
        }
      p = p->next;
    }

  return (n_letters);
}

int find_chars(const potrace_path_t * p, const Image &orig, vector<letters_t> &letters, vector<atom_t> &atom, vector<
               bond_t> &bond, int n_atom, int n_bond, int height, int width, ColorGray &bgColor, double THRESHOLD,
               int max_font_width, int max_font_height, int &real_font_width, int &real_font_height, bool verbose)
{
  int n, *tag, n_letters = 0;
  potrace_dpoint_t (*c)[3];
  real_font_width = 0;
  real_font_height = 0;

  while (p != NULL)
    {
      if ((p->sign == int('+')))
        {
          n = p->curve.n;
          tag = p->curve.tag;
          c = p->curve.c;
          int top = height;
          int x1 = 0;
          int left = width;
          int y1 = 0;
          int bottom = 0;
          int x2 = 0;
          int right = 0;
          int y2 = 0;
          int cx,cy;
          for (int i = 0; i < n; i++)
            {
              switch (tag[i])
                {
                case POTRACE_CORNER:
                  cx = c[i][1].x;
                  cy = c[i][1].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;

                  if (cx < left)
                    {
                      left = cx;
                      y1 = cy;
                    }
                  if (cx > right)
                    {
                      right = cx;
                      y2 = cy;
                    }
                  if (cy < top)
                    {
                      top = cy;
                      x1 = cx;
                    }
                  if (cy > bottom)
                    {
                      bottom = cy;
                      x2 = cx;
                    }
                  break;
                case POTRACE_CURVETO:
                  cx = c[i][0].x;
                  cy = c[i][0].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  if (cx < left)
                    {
                      left = cx;
                      y1 = cy;
                    }
                  if (cx > right)
                    {
                      right = cx;
                      y2 = cy;
                    }
                  if (cy < top)
                    {
                      top = cy;
                      x1 = cx;
                    }
                  if (cy > bottom)
                    {
                      bottom = cy;
                      x2 = cx;
                    }
                  cx = c[i][1].x;
                  cy = c[i][1].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  if (cx < left)
                    {
                      left = cx;
                      y1 = cy;
                    }
                  if (cx > right)
                    {
                      right = cx;
                      y2 = cy;
                    }
                  if (cy < top)
                    {
                      top = cy;
                      x1 = cx;
                    }
                  if (cy > bottom)
                    {
                      bottom = cy;
                      x2 = cx;
                    }
                  break;
                }
              cx = c[i][2].x;
              cy = c[i][2].y;
              if (cx<0) cx=0;
              if (cx>width) cx=width;
              if (cy<0) cy=0;
              if (cy>height) cy=height;
              if (cx < left)
                {
                  left = cx;
                  y1 = cy;
                }
              if (cx > right)
                {
                  right = cx;
                  y2 = cy;
                }
              if (cy < top)
                {
                  top = cy;
                  x1 = cx;
                }
              if (cy > bottom)
                {
                  bottom = cy;
                  x2 = cx;
                }
            }

          if (((bottom - top) <= 2 * max_font_height) && ((right - left) <= 2 * max_font_width) && (right - left
              > V_DISPLACEMENT) && (bottom - top > MIN_FONT_HEIGHT))
            {
              int s = 1;
              while ((top > 0) && (s > 0))
                {
                  s = 0;
                  s = get_pixel(orig, bgColor, x1, top, THRESHOLD);
                  if (s > 0)
                    top--;
                }
              s = 1;
              while ((bottom < height) && (s > 0))
                {
                  s = 0;
                  s = get_pixel(orig, bgColor, x2, bottom, THRESHOLD);
                  if (s > 0)
                    bottom++;
                }
              s = 1;
              while ((left > 0) && (s > 0))
                {
                  s = 0;
                  s = get_pixel(orig, bgColor, left, y1, THRESHOLD);
                  if (s > 0)
                    left--;
                }
              s = 1;
              while ((right < width) && (s > 0))
                {
                  s = 0;
                  s = get_pixel(orig, bgColor, right, y2, THRESHOLD);
                  if (s > 0)
                    right++;
                }
            }

          bool found = false;
          if (((bottom - top) <= max_font_height) && ((right - left) <= max_font_width) && (right - left
              > V_DISPLACEMENT) && (bottom - top > MIN_FONT_HEIGHT))
            {

              char label = 0;
              label = get_atom_label(orig, bgColor, left, top, right, bottom, THRESHOLD, (right + left) / 2, top, false, verbose);

              if (label != 0)
                {
                  letters_t lt;
                  letters.push_back(lt);
                  letters[n_letters].a = label;
                  letters[n_letters].x = (left + right) / 2;
                  letters[n_letters].y = (top + bottom) / 2;
                  letters[n_letters].r = distance(left, top, right, bottom) / 2;
		  letters[n_letters].min_x = left;
		  letters[n_letters].max_x = right;
		  letters[n_letters].min_y = top;
		  letters[n_letters].max_y = bottom;
		  letters[n_letters].curve = p;
                  if (right - left > real_font_width)
                    real_font_width = right - left;
                  if (bottom - top > real_font_height)
                    real_font_height = bottom - top;
                  letters[n_letters].free = true;
                  n_letters++;
                  if (n_letters >= MAX_ATOMS)
                    n_letters--;
                  delete_bonds_in_char(bond, n_bond, atom, left, top, right, bottom);
                  delete_curve_with_children(atom, bond, n_atom, n_bond, p);
                  found = true;
                }
            }
          if (((bottom - top) <= 2 * max_font_height) && ((right - left) <= max_font_width) && (right - left
              > V_DISPLACEMENT) && (bottom - top > MIN_FONT_HEIGHT) && !found)
            {

              char label1 = 0;
              int newtop = (top + bottom) / 2;
              label1 = get_atom_label(orig, bgColor, left, newtop, right, bottom, THRESHOLD, (right + left) / 2,
                                      newtop, false, verbose);
              char label2 = 0;
              int newbottom = (top + bottom) / 2;
              label2 = get_atom_label(orig, bgColor, left, top, right, newbottom, THRESHOLD, (right + left) / 2, top, false, verbose);
              if (label1 != 0 && label2 != 0 && !(tolower(label1) == 's' && tolower(label2) =='s'))
                {
                  //cout << label1 << label2 << endl;
                  letters_t lt1;
                  letters.push_back(lt1);
                  letters[n_letters].a = label1;
                  letters[n_letters].x = (left + right) / 2;
                  letters[n_letters].y = (newtop + bottom) / 2;
                  letters[n_letters].r = distance(left, newtop, right, bottom) / 2;
		  letters[n_letters].min_x = left;
		  letters[n_letters].max_x = right;
		  letters[n_letters].min_y = newtop;
		  letters[n_letters].max_y = bottom;
		  letters[n_letters].curve = p;
                  if (right - left > real_font_width)
                    real_font_width = right - left;
                  if (bottom - newtop > real_font_height)
                    real_font_height = bottom - newtop;
                  letters[n_letters].free = true;
                  n_letters++;
                  if (n_letters >= MAX_ATOMS)
                    n_letters--;
                  letters_t lt2;
                  letters.push_back(lt2);
                  letters[n_letters].a = label2;
                  letters[n_letters].x = (left + right) / 2;
                  letters[n_letters].y = (top + newbottom) / 2;
                  letters[n_letters].r = distance(left, top, right, newbottom) / 2;
		  letters[n_letters].min_x = left;
		  letters[n_letters].max_x = right;
		  letters[n_letters].min_y = top;
		  letters[n_letters].max_y = newbottom;
		  letters[n_letters].curve = p;
                  if (newbottom - top > real_font_height)
                    real_font_height = newbottom - top;
                  letters[n_letters].free = true;
                  n_letters++;
                  if (n_letters >= MAX_ATOMS)
                    n_letters--;
                  delete_bonds_in_char(bond, n_bond, atom, left, top, right, bottom);
                  delete_curve_with_children(atom, bond, n_atom, n_bond, p);
                  found = true;
                }
            }
          if (((bottom - top) <= max_font_height) && ((right - left) <= 2 * max_font_width) && (right - left
              > V_DISPLACEMENT) && (bottom - top > MIN_FONT_HEIGHT) && !found)
            {

              char label1 = 0;
              int newright = (left + right) / 2;
              label1 = get_atom_label(orig, bgColor, left, top, newright, bottom, THRESHOLD, (left + newright) / 2,
                                      top, false, verbose);
              char label2 = 0;
              int newleft = (left + right) / 2;
              label2 = get_atom_label(orig, bgColor, newleft, top, right, bottom, THRESHOLD, (newleft + right) / 2,
                                      top, false, verbose);
              if ((label1 != 0) && (label2 != 0))
                {
                  //cout << label1 << label2 << endl;
                  letters_t lt1;
                  letters.push_back(lt1);
                  letters[n_letters].a = label1;
                  letters[n_letters].x = (left + newright) / 2;
                  letters[n_letters].y = (top + bottom) / 2;
                  letters[n_letters].r = distance(left, top, newright, bottom) / 2;
		  letters[n_letters].min_x = left;
		  letters[n_letters].max_x = newright;
		  letters[n_letters].min_y = top;
		  letters[n_letters].max_y = bottom;
		  letters[n_letters].curve = p;
                  if (newright - left > real_font_width)
                    real_font_width = newright - left;
                  if (bottom - top > real_font_height)
                    real_font_height = bottom - top;
                  letters[n_letters].free = true;
                  n_letters++;
                  if (n_letters >= MAX_ATOMS)
                    n_letters--;
                  letters_t lt2;
                  letters.push_back(lt2);
                  letters[n_letters].a = label2;
                  letters[n_letters].x = (newleft + right) / 2;
                  letters[n_letters].y = (top + bottom) / 2;
                  letters[n_letters].r = distance(newleft, top, right, bottom) / 2;
		  letters[n_letters].min_x = newleft;
		  letters[n_letters].max_x = right;
		  letters[n_letters].min_y = top;
		  letters[n_letters].max_y = bottom;
		  letters[n_letters].curve = p;
                  if (right - newleft > real_font_width)
                    real_font_width = right - newleft;
                  letters[n_letters].free = true;
                  n_letters++;
                  if (n_letters >= MAX_ATOMS)
                    n_letters--;
                  delete_bonds_in_char(bond, n_bond, atom, left, top, right, bottom);
                  delete_curve_with_children(atom, bond, n_atom, n_bond, p);
                  found = true;
                }
            }

        }
      p = p->next;
    }
  if (real_font_width < 1)
    real_font_width = max_font_width;
  else
    real_font_width++;
  if (real_font_height < 1)
    real_font_height = max_font_height;
  else
    real_font_height++;
  return (n_letters);
}



int find_fused_chars(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, vector<letters_t> &letters, int n_letters,
                     int max_font_height, int max_font_width, char dummy, const Image &orig, const ColorGray &bgColor,
                     double THRESHOLD, unsigned int size, bool verbose)
{
  double dist = max(max_font_width, max_font_height);

  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && bond_length(bond, i, atom) < dist)
      {
        list<int> t;
        t.push_back(i);
        double xmin1 = min(atom[bond[i].a].x, atom[bond[i].b].x);
        double xmax1 = max(atom[bond[i].a].x, atom[bond[i].b].x);
        double ymin1 = min(atom[bond[i].a].y, atom[bond[i].b].y);
        double ymax1 = max(atom[bond[i].a].y, atom[bond[i].b].y);
        for (int j = 0; j < n_bond; j++)
          if (bond[j].exists && bond_length(bond, j, atom) < dist && j != i && atom[bond[j].a].x >= xmin1
              && atom[bond[j].a].x >= xmin1)
            {
              double xmax2 = max(xmax1, max(atom[bond[j].a].x, atom[bond[j].b].x));
              double ymin2 = min(ymin1, min(atom[bond[j].a].y, atom[bond[j].b].y));
              double ymax2 = max(ymax1, max(atom[bond[j].a].y, atom[bond[j].b].y));

              if (xmax2 - xmin1 <= max_font_width && ymax2 - ymin2 <= max_font_height)
                t.push_back(j);
            }

        vector<int> all_bonds(n_bond, 0);
        for (int j = 0; j < n_bond; j++)
          if (bond[j].exists)
            all_bonds[j] = 1;

        list<int> bag1, bag2;
        all_bonds[i] = 2;
        bag1.push_back(i);
        while (!bag1.empty())
          {
            int k = bag1.front();
            bag1.pop_front();
            all_bonds[k] = 3;
            for (int j = 0; j < n_bond; j++)
              if (j != k && all_bonds[j] == 1 && (bond[k].a == bond[j].a || bond[k].a == bond[j].b || bond[k].b
                                                  == bond[j].a || bond[k].b == bond[j].b))
                {
                  all_bonds[j] = 2;
                  bag1.push_back(j);
                }
          }
        while (!t.empty())
          {
            int k = t.front();
            t.pop_front();
            if (all_bonds[k] == 3)
              bag2.push_back(k);
          }

        unsigned int bag_size = bag2.size();
        if (bag_size > size)
          {
            double cx = 0;
            double cy = 0;
            int n = 0;
            double l = FLT_MAX, r = 0, t = FLT_MAX, b = 0;

            while (!bag2.empty())
              {
                int k = bag2.front();
                bag2.pop_front();
                cx += atom[bond[k].a].x + atom[bond[k].b].x;
                cy += atom[bond[k].a].y + atom[bond[k].b].y;
                l = min(l, min(atom[bond[k].a].x, atom[bond[k].b].x));
                r = max(r, max(atom[bond[k].a].x, atom[bond[k].b].x));
                t = min(t, min(atom[bond[k].a].y, atom[bond[k].b].y));
                b = max(b, max(atom[bond[k].a].y, atom[bond[k].b].y));
                n += 2;
              }
            cx /= n;
            cy /= n;
            if (r - l > MIN_FONT_HEIGHT && b - t > MIN_FONT_HEIGHT)
              {
                int left = int(cx - max_font_width / 2) - 1;
                int right = int(cx + max_font_width / 2) - 1;
                int top = int(cy - max_font_height / 2);
                int bottom = int(cy + max_font_height / 2);
                char label = 0;
                if (dummy != 0)
                  {

                    label = dummy;
                    //cout << bag_size << " " << left << " " << top << " " << right << " " << bottom << endl;
                  }
                else
                  {
                    label = get_atom_label(orig, bgColor, left, top, right, bottom, THRESHOLD, (left + right) / 2,
                                           top, false, verbose);
                  }
                if ((label != 0 && label != 'P' && label != 'p' && label != 'F' && label != 'X' && label != 'Y'
                     && label != 'n' && label != 'F' && label != 'U' && label != 'u' && label != 'h') || dummy
                    != 0)
                  {

                    bool overlap = false;
                    for (int j = 0; j < n_letters; j++)
                      {
                        if (distance((left + right) / 2, (top + bottom) / 2, letters[j].x, letters[j].y)
                            < letters[j].r)
                          overlap = true;
                      }
                    if (!overlap)
                      {
                        letters_t lt;
                        letters.push_back(lt);
                        letters[n_letters].a = label;
                        letters[n_letters].x = (left + right) / 2;
                        letters[n_letters].y = (top + bottom) / 2;
                        letters[n_letters].r = distance(left, top, right, bottom) / 2;
			letters[n_letters].min_x = l;
			letters[n_letters].max_x = r;
			letters[n_letters].min_y = t;
			letters[n_letters].max_y = b;
                        letters[n_letters].free = true;
                        n_letters++;
                        if (n_letters >= MAX_ATOMS)
                          n_letters--;
                      }
                    delete_bonds_in_char(bond, n_bond, atom, left, top, right, bottom);
                  }

              }
          }
      }
  return (n_letters);
}

void detect_plus_minus(const Image &image, ColorGray &bg, double THRESHOLD, int x1, int x2, int y1, int y2, int top, int left, int right, int bottom, bool &is_plus, bool &is_minus)
{
  vector < vector <short> > pic(right-left+1, vector<short> (bottom-top+1,0));
  for (int i=left; i<=right; i++)
    for (int j=top; j<=bottom; j++)
      pic[i-left][j-top] = get_pixel(image, bg, i, j, THRESHOLD) ;

  int y = 0;
  int x = x1 - left;
  if (x<0 || x>=pic.size()) return;
  while (y<pic[x].size() &&  pic[x][y] != 1) y++;
  if (y>=pic[x].size()) return;
  vector<point_t> points;
  // populating points with BFS
  list<point_t> bag;
  point_t p;
  p.x = x;
  p.y = y;
  bag.push_back(p);

  pic[p.x][p.y] = 2;
  
  while (!bag.empty())
    {
      point_t c = bag.front();
      bag.pop_front();
      points.push_back(c);
      pic[c.x][c.y] = 0;
        // this goes around 3x3 square touching the chosen pixel
      for (int i = c.x - 1; i < c.x + 2; i++)
	for (int j = c.y - 1; j < c.y + 2; j++)
	  if (i >= 0 && j >= 0 && i < pic.size() && j < pic[i].size()  && pic[i][j] == 1)
              {
		point_t t;
                t.x = i;
                t.y = j;
		bag.push_back(t);
                pic[i][j] = 2;
              }
      }
  const int len=50;
  vector<int> hist(len,0);
  int top_pos=0;
  int top_value=0;
  point_t head, tail,center;
  int min_x, min_y, max_x, max_y;
  build_hist(points,hist,len,top_pos,top_value,head,tail,center,min_x, min_y, max_x, max_y);
  double fill = (double) points.size() / ((max_x-min_x+1)*(max_y-min_y+1));
  double aspect = (double) (max_y - min_y +1) / (max_x - min_x +1);
  if (points.size() > len)
    {     
      vector<int> peaks(1,top_pos);
      vector<int> values(1,top_value);
      for (int k=1; k<len;k++)
	{
	  int pos=k+top_pos;
	  if (pos>=len) pos -= len;
	  int after=pos+1;
	  int before=pos-1;
	  if (after>=len) after -=len;
	  if (before<0) before +=len;
	  if (hist[before]<hist[pos] && hist[after]<hist[pos] && hist[pos]>=top_value/2)  // find all peaks at least half as high as the top-most
	    {
	      peaks.push_back(pos);
	      values.push_back(hist[pos]);
	    }
	}

      if (peaks.size() == 2   && abs(len/2 - abs(peaks[1]-peaks[0]))<=1)  // only two peaks are present at 180 degrees 
	{
	  
	  if (aspect < 0.7 && fill > 0.9)
	    is_minus = true;
	}
      
      if (peaks.size() == 4  && (double(values[1])/values[0]>0.8 || values[0]-values[1]<=2)  && (double(values[2])/values[0]>0.8  || values[0]-values[2]<=2) && (double(values[3])/values[0]>0.8 || values[0]-values[3]<=2))
	{
	  bool first=false, second=false, third=false, fourth=false;
	  for (int j=0; j<4; j++)
	    {
	      if (peaks[j] <= 1 || peaks[j] >=len-2) first=true;
	      if (abs(len/4-peaks[j])<=1) second=true;
	      if (abs(len/2-peaks[j])<=1) third=true;
	      if (abs(3*len/4-peaks[j])<=1) fourth=true;
	    }
	  for (int j=0; j<peaks.size(); j++)          // check outside of the peaks is essentially zero
	    for (int k=peaks[j]-2; k<=peaks[j]+2; k++)
	      {
		int kk=k;
		if (kk<0) kk += len;
		if (kk>=len) kk -=len;
		hist[kk]=0;
	      }
	  
	  bool low=true;
	  for(int k=0; k<len; k++)
	    if (hist[k]>3) low=false;
	  if (first && second && third && fourth && low)
	    {
	      // we found a plus!
	      is_plus = true;
	    }
	}	
    }
  else
    {
      if (aspect < 0.7 && fill > 0.9)
	is_minus = true;
      if (aspect > 0.7 && aspect < 1. / 0.7 && abs(y1 - y2) < 3 && abs(y1 + y2 - bottom - top) / 2 < 3
	  && abs(x1 - x2) < 3 && abs(x1 + x2 - right - left) / 2 < 3)
	is_plus = true;
    }
  
}

int find_plus_minus(const potrace_path_t *p, const Image &image, ColorGray &bgColor, double THRESHOLD, vector<letters_t> &letters, vector<atom_t> &atom, vector<bond_t> &bond,
                    int n_atom, int n_bond, int height, int width, int max_font_height, int max_font_width, int n_letters, double avg_bond_length)
{
  int n, *tag;
  potrace_dpoint_t (*c)[3];

  while (p != NULL)
    {
      if ((p->sign == int('+')) && detect_curve(bond, n_bond, p))
        {
          n = p->curve.n;
          tag = p->curve.tag;
          c = p->curve.c;
          int top = height;
          int x1 = 0;
          int left = width;
          int y1 = 0;
          int bottom = 0;
          int x2 = 0;
          int right = 0;
          int y2 = 0;
          int cx,cy;
          for (int i = 0; i < n; i++)
            {
              switch (tag[i])
                {
                case POTRACE_CORNER:
                  cx = c[i][1].x;
                  cy = c[i][1].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  if (cx < left)
                    {
                      left = cx;
                      y1 = cy;
                    }
                  if (cx > right)
                    {
                      right = cx;
                      y2 = cy;
                    }
                  if (cy < top)
                    {
                      top = cy;
                      x1 = cx;
                    }
                  if (cy > bottom)
                    {
                      bottom = cy;
                      x2 = cx;
                    }
                  break;
                case POTRACE_CURVETO:
                  cx = c[i][0].x;
                  cy = c[i][0].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  if (cx < left)
                    {
                      left = cx;
                      y1 = cy;
                    }
                  if (cx > right)
                    {
                      right = cx;
                      y2 = cy;
                    }
                  if (cy < top)
                    {
                      top = cy;
                      x1 = cx;
                    }
                  if (cy > bottom)
                    {
                      bottom = cy;
                      x2 = cx;
                    }
                  cx = c[i][1].x;
                  cy = c[i][1].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  if (cx < left)
                    {
                      left = cx;
                      y1 = cy;
                    }
                  if (cx > right)
                    {
                      right = cx;
                      y2 = cy;
                    }
                  if (cy < top)
                    {
                      top = cy;
                      x1 = cx;
                    }
                  if (cy > bottom)
                    {
                      bottom = cy;
                      x2 = cx;
                    }
                  break;
                }
              cx = c[i][2].x;
              cy = c[i][2].y;
              if (cx<0) cx=0;
              if (cx>width) cx=width;
              if (cy<0) cy=0;
              if (cy>height) cy=height;
              if (cx < left)
                {
                  left = cx;
                  y1 = cy;
                }
              if (cx > right)
                {
                  right = cx;
                  y2 = cy;
                }
              if (cy < top)
                {
                  top = cy;
                  x1 = cx;
                }
              if (cy > bottom)
                {
                  bottom = cy;
                  x2 = cx;
                }
            }

          if (((bottom - top) <= max_font_height) && ((right - left) <= max_font_width) && (right - left > 1))
            {
              char c = ' ';
              bool char_to_right = false;
              bool inside_char = false;
              for (int j = 0; j < n_letters; j++)
                {
                  if (letters[j].x > right && (top + bottom) / 2 > letters[j].min_y  && (top + bottom) / 2 < letters[j].max_y 
		      && right > letters[j].min_x - letters[j].r && letters[j].a  != '-' && letters[j].a != '+')
		    char_to_right = true;
                  if (letters[j].min_x <= left && letters[j].max_x >= right && letters[j].min_y <= top && letters[j].max_y >= bottom)
                    inside_char = true;
                }
              bool is_minus = false;
	      bool is_plus = false;
	      detect_plus_minus(image,bgColor, THRESHOLD,x1,x2,y1,y2,top,left,right,bottom,is_plus,is_minus);

              if (is_minus && !char_to_right && !inside_char && (right - left) < avg_bond_length/2)
		{
		  c = '-';
		  //		  cout<<"Minus: "<<(right-left)<<" "<<max_font_width<<endl;
		}
              if (is_plus && !inside_char)
		c = '+';

              if (c != ' ')
                {
                  letters_t lt;
                  letters.push_back(lt);
                  letters[n_letters].a = c;
                  letters[n_letters].x = (left + right) / 2;
                  letters[n_letters].y = (top + bottom) / 2;
                  letters[n_letters].r = distance(left, top, right, bottom) / 2;
		  letters[n_letters].min_x = left;
		  letters[n_letters].max_x = right;
		  letters[n_letters].min_y = top;
		  letters[n_letters].max_y = bottom;
                  letters[n_letters].free = true;
                  n_letters++;
                  if (n_letters >= MAX_ATOMS)
                    n_letters--;
                  delete_curve_with_children(atom, bond, n_atom, n_bond, p);
                }
            }
        }
      p = p->next;
    }
  return (n_letters);
}

int clean_unrecognized_characters(vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom, int real_font_height,
                                  int real_font_width, unsigned int size, vector<letters_t> &letters, int n_letters)
{
  vector<int> all_bonds(n_bond, 0);

  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      all_bonds[i] = 1;

  for (int i = 0; i < n_bond; i++)
    if (all_bonds[i] == 1)
      {
        list<int> bag1, bag2;
        list<int> trash;
        all_bonds[i] = 2;
        bag1.push_back(i);
        while (!bag1.empty())
          {
            int k = bag1.front();
            bag1.pop_front();
            all_bonds[k] = 0;
            bag2.push_back(k);
            for (int j = 0; j < n_bond; j++)
              if (j != k && all_bonds[j] == 1 && (bond[k].a == bond[j].a || bond[k].a == bond[j].b || bond[k].b
                                                  == bond[j].a || bond[k].b == bond[j].b))
                {
                  all_bonds[j] = 2;
                  bag1.push_back(j);
                }
          }
        double t = FLT_MAX, b = 0, l = FLT_MAX, r = 0;
        while (!bag2.empty())
          {
            int k = bag2.front();
            bag2.pop_front();
            trash.push_back(k);

            if (atom[bond[k].a].x < l)
              l = atom[bond[k].a].x;
            if (atom[bond[k].b].x < l)
              l = atom[bond[k].b].x;
            if (atom[bond[k].a].x > r)
              r = atom[bond[k].a].x;
            if (atom[bond[k].b].x > r)
              r = atom[bond[k].b].x;
            if (atom[bond[k].a].y < t)
              t = atom[bond[k].a].y;
            if (atom[bond[k].b].y < t)
              t = atom[bond[k].b].y;
            if (atom[bond[k].a].y > b)
              b = atom[bond[k].a].y;
            if (atom[bond[k].b].y > b)
              b = atom[bond[k].b].y;
          }
        if ((r - l) < real_font_width && (b - t) < real_font_height && trash.size() > size)
          {
            while (!trash.empty())
              {
                int k = trash.front();
                trash.pop_front();
                bond[k].exists = false;
              }
            letters_t lt;
            letters.push_back(lt);
            letters[n_letters].a = '*';
            letters[n_letters].x = (l + r) / 2;
            letters[n_letters].y = (t + b) / 2;
            letters[n_letters].r = distance(l, t, r, b) / 2;
	    letters[n_letters].min_x = l;
	    letters[n_letters].max_x = r;
	    letters[n_letters].min_y = t;
	    letters[n_letters].max_y = b;
            letters[n_letters].free = true;
	    // cout<<n_letters<<" "<<l<<" "<<t<<" "<<r<<" "<<b<<endl;
            n_letters++;
            if (n_letters >= MAX_ATOMS)
              n_letters--;
          }
      }
  return (n_letters);
}

void remove_small_terminal_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double avg)
{
  bool found = true;

  while (found)
    {
      found = false;
      for (int j = 0; j < n_bond; j++)
        if (bond[j].exists && bond[j].type == 1 && !bond[j].wedge && !bond[j].hash && !bond[j].arom && bond_length(bond, j, atom) < avg / 2)
          {
            bool not_corner_a = terminal_bond(bond[j].a, j, bond, n_bond);
            bool not_corner_b = terminal_bond(bond[j].b, j, bond, n_bond);
            if (not_corner_a)
              {
                bond[j].exists = false;
                atom[bond[j].a].exists = false;
                found = true;
                if (atom[bond[j].b].label == " ")
                  {
                    if (atom[bond[j].a].label != " ")
		      {
			atom[bond[j].b].label = atom[bond[j].a].label;
			atom[bond[j].b].min_x = min(atom[bond[j].b].min_x,atom[bond[j].a].min_x);
			atom[bond[j].b].min_y = min(atom[bond[j].b].min_y,atom[bond[j].a].min_y);
			atom[bond[j].b].max_x = max(atom[bond[j].b].max_x,atom[bond[j].a].max_x);
			atom[bond[j].b].max_y = max(atom[bond[j].b].max_y,atom[bond[j].a].max_y);
		      }
                    else
                      {
                        bool dashed = false;
                        int n = 0;
                        for (int i = 0; i < n_bond; i++)
                          if (bond[i].exists && i != j && (bond[i].a == bond[j].b || bond[i].b == bond[j].b))
                            {
                              n++;
                              if (bond[i].hash)
                                dashed = true;
                            }
                        if (!dashed)
			  {
			    atom[bond[j].b].label = "Xx";
			    atom[bond[j].b].min_x = min(atom[bond[j].b].min_x,atom[bond[j].a].min_x);
			    atom[bond[j].b].min_y = min(atom[bond[j].b].min_y,atom[bond[j].a].min_y);
			    atom[bond[j].b].max_x = max(atom[bond[j].b].max_x,atom[bond[j].a].max_x);
			    atom[bond[j].b].max_y = max(atom[bond[j].b].max_y,atom[bond[j].a].max_y);
			  }
                      }
                  }
              }
            if (not_corner_b)
              {
                bond[j].exists = false;
                atom[bond[j].b].exists = false;
                found = true;
                if (atom[bond[j].a].label == " ")
                  {
                    if (atom[bond[j].b].label != " ")
		      {
			atom[bond[j].a].label = atom[bond[j].b].label;
			atom[bond[j].a].min_x = min(atom[bond[j].b].min_x,atom[bond[j].a].min_x);
			atom[bond[j].a].min_y = min(atom[bond[j].b].min_y,atom[bond[j].a].min_y);
			atom[bond[j].a].max_x = max(atom[bond[j].b].max_x,atom[bond[j].a].max_x);
			atom[bond[j].a].max_y = max(atom[bond[j].b].max_y,atom[bond[j].a].max_y);
		      }
                    else
                      {
                        bool dashed = false;
                        int n = 0;
                        for (int i = 0; i < n_bond; i++)
                          if (bond[i].exists && i != j && (bond[i].a == bond[j].a || bond[i].b == bond[j].a))
                            {
                              n++;
                              if (bond[i].hash)
                                dashed = true;
                            }
                        if (!dashed)
			  {
			    atom[bond[j].a].label = "Xx";
			    atom[bond[j].a].min_x = min(atom[bond[j].b].min_x,atom[bond[j].a].min_x);
			    atom[bond[j].a].min_y = min(atom[bond[j].b].min_y,atom[bond[j].a].min_y);
			    atom[bond[j].a].max_x = max(atom[bond[j].b].max_x,atom[bond[j].a].max_x);
			    atom[bond[j].a].max_y = max(atom[bond[j].b].max_y,atom[bond[j].a].max_y);
			  }
                      }
                  }
              }
          }
    }
}
