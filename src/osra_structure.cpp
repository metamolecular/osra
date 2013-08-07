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

// File: osra_structure.cpp
//
// Defines main structure recognition (molecular atoms and bonds)
// recognition routines
//

#include <math.h> // fabs(double)
#include <float.h> // FLT_MAX
#include <iostream> // std::cout
#include <set>

#include "osra_common.h"
#include "osra_structure.h"
#include "osra_ocr.h"
#include "osra_openbabel.h"

void remove_disconnected_atoms(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond)
{
  for (int i = 0; i < n_atom; i++)
    {
      if (atom[i].exists)
        {
          atom[i].exists = false;
          for (int j = 0; j < n_bond; j++)
            {
              if ((bond[j].exists) && (i == bond[j].a || i == bond[j].b))
                {
                  atom[i].exists = true;
                }
            }
        }
    }
}

void remove_zero_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom)
{
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      {
        for (int j = 0; j < n_bond; j++)
          if ((bond[j].exists) && (j != i) && ((bond[i].a == bond[j].a && bond[i].b == bond[j].b) || (bond[i].a == bond[j].b && bond[i].b == bond[j].a)))
            bond[j].exists = false;
        if (bond[i].a == bond[i].b)
          bond[i].exists = false;
        if (!atom[bond[i].a].exists || !atom[bond[i].b].exists)
          bond[i].exists = false;
      }
}

void collapse_doubleup_bonds(vector<bond_t> &bond, int n_bond)
{
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      for (int j = 0; j < n_bond; j++)
        if ((bond[j].exists) && (j != i) && ((bond[i].a == bond[j].a && bond[i].b == bond[j].b) || (bond[i].a == bond[j].b && bond[i].b == bond[j].a)))
          {
            bond[j].exists = false;
            bond[i].type++;
          }
}

void bond_end_swap(vector<bond_t> &bond, int i)
{
  int t = bond[i].a;
  bond[i].a = bond[i].b;
  bond[i].b = t;
}

bool bonds_within_each_other(const vector<bond_t> &bond, int ii, int jj, const vector<atom_t> &atom)
{
  int i, j;
  bool res = false;

  if (bond_length(bond, ii, atom) > bond_length(bond, jj, atom))
    {
      i = ii;
      j = jj;
    }
  else
    {
      i = jj;
      j = ii;
    }

  double x1 = atom[bond[i].a].x;
  double x2 = atom[bond[i].b].x;
  double y1 = atom[bond[i].a].y;
  double y2 = atom[bond[i].b].y;
  double d1 = bond_length(bond, i, atom);
  double x3 = distance_from_bond_x_a(x1, y1, x2, y2, atom[bond[j].a].x, atom[bond[j].a].y);
  double x4 = distance_from_bond_x_a(x1, y1, x2, y2, atom[bond[j].b].x, atom[bond[j].b].y);

  if ((x3 + x4) / 2 > 0 && (x3 + x4) / 2 < d1)
    res = true;

  return (res);
}

bool no_white_space(int ai, int bi, int aj, int bj, const vector<atom_t> &atom, const Image &image, double threshold,
                    const ColorGray &bgColor)
{
  vector<double> xx(4);
  double dx1 = atom[bi].x - atom[ai].x;
  double dy1 = atom[bi].y - atom[ai].y;
  double dx2 = atom[bj].x - atom[aj].x;
  double dy2 = atom[bj].y - atom[aj].y;
  double k1, k2;
  int s = 0, w = 0;
  int total_length = 0, white_length = 0;

  if (fabs(dx1) > fabs(dy1))
    {
      xx[0] = atom[ai].x;
      xx[1] = atom[bi].x;
      xx[2] = atom[aj].x;
      xx[3] = atom[bj].x;
      std::sort(xx.begin(), xx.end());
      k1 = dy1 / dx1;
      k2 = dy2 / dx2;
      int d = (dx1 > 0 ? 1 : -1);

      for (int x = int(atom[ai].x); x != int(atom[bi].x); x += d)
        if (x > xx[1] && x < xx[2])
          {
            double p1 = (x - atom[ai].x) * k1 + atom[ai].y;
            double p2 = (x - atom[aj].x) * k2 + atom[aj].y;
            if (fabs(p2 - p1) < 1)
              continue;
            int dp = (p2 > p1 ? 1 : -1);
            bool white = false;
            for (int y = int(p1) + dp; y != int(p2); y += dp)
              {
                s++;
                if (get_pixel(image, bgColor, x, y, threshold) == 0)
                  {
                    w++;
                    white = true;
                  }
              }
            total_length++;
            if (white)
              white_length++;
          }
    }
  else
    {
      xx[0] = atom[ai].y;
      xx[1] = atom[bi].y;
      xx[2] = atom[aj].y;
      xx[3] = atom[bj].y;
      std::sort(xx.begin(), xx.end());
      k1 = dx1 / dy1;
      k2 = dx2 / dy2;
      int d = (dy1 > 0 ? 1 : -1);

      for (int y = int(atom[ai].y); y != int(atom[bi].y); y += d)
        if (y > xx[1] && y < xx[2])
          {
            double p1 = (y - atom[ai].y) * k1 + atom[ai].x;
            double p2 = (y - atom[aj].y) * k2 + atom[aj].x;
            if (fabs(p2 - p1) < 1)
              continue;
            int dp = (p2 > p1 ? 1 : -1);
            bool white = false;
            for (int x = int(p1) + dp; x != int(p2); x += dp)
              {
                s++;
                if (get_pixel(image, bgColor, x, y, threshold) == 0)
                  {
                    w++;
                    white = true;
                  }
              }
            total_length++;
            if (white)
              white_length++;
          }
    }
  //if (s == 0) return(true);
  //if ((1. * w) / s > WHITE_SPACE_FRACTION) return(false);
  if (total_length == 0)
    return (true);
  if ((1. * white_length) / total_length > 0.5)
    return (false);
  else
    return (true);

}

int find_wavy_bonds(vector<bond_t> &bond,int n_bond,const vector<atom_t> &atom,double avg)
{
  int n_bond_orig = n_bond;
  const double t = 4;
  for (int i = 0; i < n_bond_orig; i++)
    if (bond[i].exists)
      {
	double l = bond_length(bond, i, atom);
	if (l > avg)
	  avg = l;
      }

      for (int j = 0; j < n_bond_orig; j++)
	{
	  bool found = false;
	  set<int> bag;
	  int begin = 0;
	  int end = 0;
	  if (bond[j].exists)// && !bond[j].Small)
	    {
	      double l = bond_length(bond, j, atom);
	      if (l < avg/4)
		{
		  bag.insert(j);
		  begin = j;
		  end = j;
		  found = true;
		}
	    }
	  while (found)
	    {
	      found = false;
	      int mini = -1;
	      double mind = t;
	      double d1,d2,d3,d4;
	      for (int i = 0; i < n_bond_orig; i++)
		if (bond[i].exists /*&& !bond[i].Small*/ && bag.find(i)==bag.end() && 
		    ((d1 = distance(atom[bond[end].a].x,atom[bond[end].a].y,atom[bond[i].a].x,atom[bond[i].a].y)) < mind ||
		     (d2 = distance(atom[bond[end].a].x,atom[bond[end].a].y,atom[bond[i].b].x,atom[bond[i].b].y)) < mind ||
		     (d3 = distance(atom[bond[end].b].x,atom[bond[end].b].y,atom[bond[i].a].x,atom[bond[i].a].y)) < mind ||
		     (d4 = distance(atom[bond[end].b].x,atom[bond[end].b].y,atom[bond[i].b].x,atom[bond[i].b].y)) < mind ))
		  //		    (bond[end].a == bond[i].a || bond[end].a == bond[i].b || bond[end].b == bond[i].a || bond[end].b == bond[i].b))
		  {
		    mind = min(min(d1,d2),min(d3,d4));
		    mini = i;
		  }
	      if (mini >= 0)
		{
		  double l = bond_length(bond, mini, atom);
		  if (l < avg/4)
		    {

		      bag.insert(mini);
		      end = mini;
		      found = true;
		    }
		}
	      mini = -1;
	      mind = t;
	      for (int i = 0; i < n_bond_orig; i++)
		if (bond[i].exists /*&& !bond[i].Small*/ && bag.find(i)==bag.end() && 
		    ((d1 = distance(atom[bond[begin].a].x,atom[bond[begin].a].y,atom[bond[i].a].x,atom[bond[i].a].y)) < mind ||
		     (d2 = distance(atom[bond[begin].a].x,atom[bond[begin].a].y,atom[bond[i].b].x,atom[bond[i].b].y)) < mind ||
		     (d3 = distance(atom[bond[begin].b].x,atom[bond[begin].b].y,atom[bond[i].a].x,atom[bond[i].a].y)) < mind ||
		     (d4 = distance(atom[bond[begin].b].x,atom[bond[begin].b].y,atom[bond[i].b].x,atom[bond[i].b].y)) < mind ))
		  //		    (bond[begin].a == bond[i].a || bond[begin].a == bond[i].b || bond[begin].b == bond[i].a || bond[begin].b == bond[i].b))
		  {
		    mind = min(min(d1,d2),min(d3,d4));
		    mini = i;
		  }
	      if (mini >= 0)
		  {
		    double l = bond_length(bond, mini, atom);
		    if (l < avg/4)
		      {

			bag.insert(mini);
			begin = mini;
			found = true;
		      }
		  }
	    }

	  if (bag.size() > 3)
	    {
	      double d1 = distance(atom[bond[begin].a].x,atom[bond[begin].a].y,atom[bond[end].a].x,atom[bond[end].a].y);
	      double d2 = distance(atom[bond[begin].a].x,atom[bond[begin].a].y,atom[bond[end].b].x,atom[bond[end].b].y);
	      double d3 = distance(atom[bond[begin].b].x,atom[bond[begin].b].y,atom[bond[end].a].x,atom[bond[end].a].y);
	      double d4 = distance(atom[bond[begin].b].x,atom[bond[begin].b].y,atom[bond[end].b].x,atom[bond[end].b].y);
	      double d = max(max(d1,d2),max(d3,d4));

	      if (d>0.75*avg && d<1.5*avg)
		{
		  double xa = atom[bond[begin].a].x;
		  double ya = atom[bond[begin].a].y;
		  double xb = atom[bond[end].a].x;
		  double yb = atom[bond[end].a].y;
		  int a = bond[begin].a;
		  int b = bond[end].a;
		  if (d2 > d1)
		    {
		      xb = atom[bond[end].b].x;
		      yb = atom[bond[end].b].y;
		      b = bond[end].b;
		    }
		  if (d3 > max(d1,d2))
		    {
		      xa = atom[bond[begin].b].x;
		      ya = atom[bond[begin].b].y;
		      xb = atom[bond[end].a].x;
		      yb = atom[bond[end].a].y;
		      a = bond[begin].b;
		      b = bond[end].a;
		    }
		  if (d4 > max(d3,max(d1,d2)))
		    {
		      xa = atom[bond[begin].b].x;
		      ya = atom[bond[begin].b].y;
		      xb = atom[bond[end].b].x;
		      yb = atom[bond[end].b].y;
		      a = bond[begin].b;
		      b = bond[end].b;
		    }
		  bool too_far = false;
		  for(set<int>::iterator i=bag.begin(); i!=bag.end(); i++)
		    {
		      double da =  fabs(distance_from_bond_y(xa,ya,xb,yb,atom[bond[*i].a].x, atom[bond[*i].a].y));
		      double db =  fabs(distance_from_bond_y(xa,ya,xb,yb,atom[bond[*i].b].x, atom[bond[*i].b].y));

		      if (da > avg/4 || db > avg/4)
			{
			  too_far = true;
			  break;
			}
		    }
		  if (!too_far)  // This looks like a wavy bond
		    {
		      bond_t nb(a,b,bond[begin].curve);
		      nb.hash = true;                        
		      nb.wedge = true;
		      if (n_bond < MAX_ATOMS)
			{
			  bond.push_back(nb);
			  n_bond++;
			}
		      for(set<int>::iterator i=bag.begin(); i!=bag.end(); i++)
			bond[*i].exists = false;
		    }
		}
	    }
	}
  return n_bond;
}

double skeletize(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, const Image &image, double threshold,
                 const ColorGray &bgColor, double dist, double avg)
{
  double thickness = 0;
  vector<double> a;
  int n = 0;

  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && !bond[i].Small)
      {
        double l1 = bond_length(bond, i, atom);
        for (int j = 0; j < n_bond; j++)
          if (i != j && bond[j].exists && bonds_within_each_other(bond, i, j, atom) && !bond[j].Small)
            {
              double tt = distance_between_bonds(bond, i, j, atom);
              double tang = angle_between_bonds(bond, i, j, atom);
              if ((fabs(tang) > D_T_TOLERANCE && no_white_space(bond[i].a, bond[i].b, bond[j].a, bond[j].b, atom,
                   image, threshold, bgColor) && tt < MAX_BOND_THICKNESS) || tt < dist)
                {
                  double l2 = bond_length(bond, j, atom);
                  a.push_back(tt);
                  n++;
                  if (l1 < l2)
                    {
                      bond[i].exists = false;
                      bond[j].type = 1;
                      if (bond[i].arom)
                        bond[j].arom = true;
                      if (l1 > avg / 2)
                        {
                          double ay = fabs(distance_from_bond_y(atom[bond[j].a].x, atom[bond[j].a].y,
                                                                atom[bond[j].b].x, atom[bond[j].b].y, atom[bond[i].a].x, atom[bond[i].a].y));
                          double axa = fabs(distance_from_bond_x_a(atom[bond[j].a].x, atom[bond[j].a].y,
                                            atom[bond[j].b].x, atom[bond[j].b].y, atom[bond[i].a].x, atom[bond[i].a].y));
                          double axb = fabs(distance_from_bond_x_b(atom[bond[j].a].x, atom[bond[j].a].y,
                                            atom[bond[j].b].x, atom[bond[j].b].y, atom[bond[i].a].x, atom[bond[i].a].y));

                          if (tang > 0 && ay > axa)
                            {
                              atom[bond[i].a].x = (atom[bond[i].a].x + atom[bond[j].a].x) / 2;
                              atom[bond[i].a].y = (atom[bond[i].a].y + atom[bond[j].a].y) / 2;
                              atom[bond[j].a].x = (atom[bond[i].a].x + atom[bond[j].a].x) / 2;
                              atom[bond[j].a].y = (atom[bond[i].a].y + atom[bond[j].a].y) / 2;

			      atom[bond[i].a].min_x = min(atom[bond[i].a].min_x,atom[bond[j].a].min_x) ;
			      atom[bond[i].a].min_y = min(atom[bond[i].a].min_y,atom[bond[j].a].min_y) ;
			      atom[bond[i].a].max_x = min(atom[bond[i].a].max_x,atom[bond[j].a].max_x) ;
			      atom[bond[i].a].max_y = min(atom[bond[i].a].max_y,atom[bond[j].a].max_y) ;
			      atom[bond[j].a].min_x = min(atom[bond[i].a].min_x,atom[bond[j].a].min_x) ;
			      atom[bond[j].a].min_y = min(atom[bond[i].a].min_y,atom[bond[j].a].min_y) ;
			      atom[bond[j].a].max_x = min(atom[bond[i].a].max_x,atom[bond[j].a].max_x) ;
			      atom[bond[j].a].max_y = min(atom[bond[i].a].max_y,atom[bond[j].a].max_y) ;
                            }
                          if (tang < 0 && ay > axb)
                            {
                              atom[bond[i].a].x = (atom[bond[i].a].x + atom[bond[j].b].x) / 2;
                              atom[bond[i].a].y = (atom[bond[i].a].y + atom[bond[j].b].y) / 2;
                              atom[bond[j].b].x = (atom[bond[i].a].x + atom[bond[j].b].x) / 2;
                              atom[bond[j].b].y = (atom[bond[i].a].y + atom[bond[j].b].y) / 2;

			      atom[bond[i].a].min_x = min(atom[bond[i].a].min_x,atom[bond[j].b].min_x) ;
			      atom[bond[i].a].min_y = min(atom[bond[i].a].min_y,atom[bond[j].b].min_y) ;
			      atom[bond[i].a].max_x = min(atom[bond[i].a].max_x,atom[bond[j].b].max_x) ;
			      atom[bond[i].a].max_y = min(atom[bond[i].a].max_y,atom[bond[j].b].max_y) ;
			      atom[bond[j].b].min_x = min(atom[bond[i].a].min_x,atom[bond[j].b].min_x) ;
			      atom[bond[j].b].min_y = min(atom[bond[i].a].min_y,atom[bond[j].b].min_y) ;
			      atom[bond[j].b].max_x = min(atom[bond[i].a].max_x,atom[bond[j].b].max_x) ;
			      atom[bond[j].b].max_y = min(atom[bond[i].a].max_y,atom[bond[j].b].max_y) ;
                            }
                          double by = fabs(distance_from_bond_y(atom[bond[j].a].x, atom[bond[j].a].y,
                                                                atom[bond[j].b].x, atom[bond[j].b].y, atom[bond[i].b].x, atom[bond[i].b].y));
                          double bxa = fabs(distance_from_bond_x_a(atom[bond[j].a].x, atom[bond[j].a].y,
                                            atom[bond[j].b].x, atom[bond[j].b].y, atom[bond[i].b].x, atom[bond[i].b].y));
                          double bxb = fabs(distance_from_bond_x_b(atom[bond[j].a].x, atom[bond[j].a].y,
                                            atom[bond[j].b].x, atom[bond[j].b].y, atom[bond[i].b].x, atom[bond[i].b].y));

                          if (tang > 0 && by > bxb)
                            {
                              atom[bond[i].b].x = (atom[bond[i].b].x + atom[bond[j].b].x) / 2;
                              atom[bond[i].b].y = (atom[bond[i].b].y + atom[bond[j].b].y) / 2;
                              atom[bond[j].b].x = (atom[bond[i].b].x + atom[bond[j].b].x) / 2;
                              atom[bond[j].b].y = (atom[bond[i].b].y + atom[bond[j].b].y) / 2;

			      atom[bond[i].b].min_x = min(atom[bond[i].b].min_x,atom[bond[j].b].min_x) ;
			      atom[bond[i].b].min_y = min(atom[bond[i].b].min_y,atom[bond[j].b].min_y) ;
			      atom[bond[i].b].max_x = min(atom[bond[i].b].max_x,atom[bond[j].b].max_x) ;
			      atom[bond[i].b].max_y = min(atom[bond[i].b].max_y,atom[bond[j].b].max_y) ;
			      atom[bond[j].b].min_x = min(atom[bond[i].b].min_x,atom[bond[j].b].min_x) ;
			      atom[bond[j].b].min_y = min(atom[bond[i].b].min_y,atom[bond[j].b].min_y) ;
			      atom[bond[j].b].max_x = min(atom[bond[i].b].max_x,atom[bond[j].b].max_x) ;
			      atom[bond[j].b].max_y = min(atom[bond[i].b].max_y,atom[bond[j].b].max_y) ;
                            }
                          if (tang < 0 && by > bxa)
                            {
                              atom[bond[i].b].x = (atom[bond[i].b].x + atom[bond[j].a].x) / 2;
                              atom[bond[i].b].y = (atom[bond[i].b].y + atom[bond[j].a].y) / 2;
                              atom[bond[j].a].x = (atom[bond[i].b].x + atom[bond[j].a].x) / 2;
                              atom[bond[j].a].y = (atom[bond[i].b].y + atom[bond[j].a].y) / 2;

			      atom[bond[i].b].min_x = min(atom[bond[i].b].min_x,atom[bond[j].a].min_x) ;
			      atom[bond[i].b].min_y = min(atom[bond[i].b].min_y,atom[bond[j].a].min_y) ;
			      atom[bond[i].b].max_x = min(atom[bond[i].b].max_x,atom[bond[j].a].max_x) ;
			      atom[bond[i].b].max_y = min(atom[bond[i].b].max_y,atom[bond[j].a].max_y) ;
			      atom[bond[j].a].min_x = min(atom[bond[i].b].min_x,atom[bond[j].a].min_x) ;
			      atom[bond[j].a].min_y = min(atom[bond[i].b].min_y,atom[bond[j].a].min_y) ;
			      atom[bond[j].a].max_x = min(atom[bond[i].b].max_x,atom[bond[j].a].max_x) ;
			      atom[bond[j].a].max_y = min(atom[bond[i].b].max_y,atom[bond[j].a].max_y) ;
                            }
                        }
                      break;
                    }
                  else
                    {
                      bond[j].exists = false;
                      bond[i].type = 1;
                      if (bond[j].arom)
                        bond[i].arom = true;
                      if (l2 > avg / 2)
                        {
                          double ay = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y,
                                                                atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[j].a].x, atom[bond[j].a].y));
                          double axa = fabs(distance_from_bond_x_a(atom[bond[i].a].x, atom[bond[i].a].y,
                                            atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[j].a].x, atom[bond[j].a].y));
                          double axb = fabs(distance_from_bond_x_b(atom[bond[i].a].x, atom[bond[i].a].y,
                                            atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[j].a].x, atom[bond[j].a].y));

                          if (tang > 0 && ay > axa)
                            {
                              atom[bond[i].a].x = (atom[bond[i].a].x + atom[bond[j].a].x) / 2;
                              atom[bond[i].a].y = (atom[bond[i].a].y + atom[bond[j].a].y) / 2;
                              atom[bond[j].a].x = (atom[bond[i].a].x + atom[bond[j].a].x) / 2;
                              atom[bond[j].a].y = (atom[bond[i].a].y + atom[bond[j].a].y) / 2;

			      atom[bond[i].a].min_x = min(atom[bond[i].a].min_x,atom[bond[j].a].min_x) ;
			      atom[bond[i].a].min_y = min(atom[bond[i].a].min_y,atom[bond[j].a].min_y) ;
			      atom[bond[i].a].max_x = min(atom[bond[i].a].max_x,atom[bond[j].a].max_x) ;
			      atom[bond[i].a].max_y = min(atom[bond[i].a].max_y,atom[bond[j].a].max_y) ;
			      atom[bond[j].a].min_x = min(atom[bond[i].a].min_x,atom[bond[j].a].min_x) ;
			      atom[bond[j].a].min_y = min(atom[bond[i].a].min_y,atom[bond[j].a].min_y) ;
			      atom[bond[j].a].max_x = min(atom[bond[i].a].max_x,atom[bond[j].a].max_x) ;
			      atom[bond[j].a].max_y = min(atom[bond[i].a].max_y,atom[bond[j].a].max_y) ;
                            }
                          if (tang < 0 && ay > axb)
                            {
                              atom[bond[j].a].x = (atom[bond[j].a].x + atom[bond[i].b].x) / 2;
                              atom[bond[j].a].y = (atom[bond[j].a].y + atom[bond[i].b].y) / 2;
                              atom[bond[i].b].x = (atom[bond[j].a].x + atom[bond[i].b].x) / 2;
                              atom[bond[i].b].y = (atom[bond[j].a].y + atom[bond[i].b].y) / 2;

			      atom[bond[i].a].min_x = min(atom[bond[i].a].min_x,atom[bond[j].b].min_x) ;
			      atom[bond[i].a].min_y = min(atom[bond[i].a].min_y,atom[bond[j].b].min_y) ;
			      atom[bond[i].a].max_x = min(atom[bond[i].a].max_x,atom[bond[j].b].max_x) ;
			      atom[bond[i].a].max_y = min(atom[bond[i].a].max_y,atom[bond[j].b].max_y) ;
			      atom[bond[j].b].min_x = min(atom[bond[i].a].min_x,atom[bond[j].b].min_x) ;
			      atom[bond[j].b].min_y = min(atom[bond[i].a].min_y,atom[bond[j].b].min_y) ;
			      atom[bond[j].b].max_x = min(atom[bond[i].a].max_x,atom[bond[j].b].max_x) ;
			      atom[bond[j].b].max_y = min(atom[bond[i].a].max_y,atom[bond[j].b].max_y) ;
                            }
                          double by = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y,
                                                                atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[j].b].x, atom[bond[j].b].y));
                          double bxa = fabs(distance_from_bond_x_a(atom[bond[i].a].x, atom[bond[i].a].y,
                                            atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[j].b].x, atom[bond[j].b].y));
                          double bxb = fabs(distance_from_bond_x_b(atom[bond[i].a].x, atom[bond[i].a].y,
                                            atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[j].b].x, atom[bond[j].b].y));

                          if (tang > 0 && by > bxb)
                            {
                              atom[bond[i].b].x = (atom[bond[i].b].x + atom[bond[j].b].x) / 2;
                              atom[bond[i].b].y = (atom[bond[i].b].y + atom[bond[j].b].y) / 2;
                              atom[bond[j].b].x = (atom[bond[i].b].x + atom[bond[j].b].x) / 2;
                              atom[bond[j].b].y = (atom[bond[i].b].y + atom[bond[j].b].y) / 2;

			      atom[bond[i].b].min_x = min(atom[bond[i].b].min_x,atom[bond[j].b].min_x) ;
			      atom[bond[i].b].min_y = min(atom[bond[i].b].min_y,atom[bond[j].b].min_y) ;
			      atom[bond[i].b].max_x = min(atom[bond[i].b].max_x,atom[bond[j].b].max_x) ;
			      atom[bond[i].b].max_y = min(atom[bond[i].b].max_y,atom[bond[j].b].max_y) ;
			      atom[bond[j].b].min_x = min(atom[bond[i].b].min_x,atom[bond[j].b].min_x) ;
			      atom[bond[j].b].min_y = min(atom[bond[i].b].min_y,atom[bond[j].b].min_y) ;
			      atom[bond[j].b].max_x = min(atom[bond[i].b].max_x,atom[bond[j].b].max_x) ;
			      atom[bond[j].b].max_y = min(atom[bond[i].b].max_y,atom[bond[j].b].max_y) ;
                            }
                          if (tang < 0 && by > bxa)
                            {
                              atom[bond[j].b].x = (atom[bond[j].b].x + atom[bond[i].a].x) / 2;
                              atom[bond[j].b].y = (atom[bond[j].b].y + atom[bond[i].a].y) / 2;
                              atom[bond[i].a].x = (atom[bond[j].b].x + atom[bond[i].a].x) / 2;
                              atom[bond[i].a].y = (atom[bond[j].b].y + atom[bond[i].a].y) / 2;

			      atom[bond[i].b].min_x = min(atom[bond[i].b].min_x,atom[bond[j].a].min_x) ;
			      atom[bond[i].b].min_y = min(atom[bond[i].b].min_y,atom[bond[j].a].min_y) ;
			      atom[bond[i].b].max_x = min(atom[bond[i].b].max_x,atom[bond[j].a].max_x) ;
			      atom[bond[i].b].max_y = min(atom[bond[i].b].max_y,atom[bond[j].a].max_y) ;
			      atom[bond[j].a].min_x = min(atom[bond[i].b].min_x,atom[bond[j].a].min_x) ;
			      atom[bond[j].a].min_y = min(atom[bond[i].b].min_y,atom[bond[j].a].min_y) ;
			      atom[bond[j].a].max_x = min(atom[bond[i].b].max_x,atom[bond[j].a].max_x) ;
			      atom[bond[j].a].max_y = min(atom[bond[i].b].max_y,atom[bond[j].a].max_y) ;
                            }
                        }
                    }
                }
            }
      }
  std::sort(a.begin(), a.end());
  if (n > 0)
    thickness = a[(n - 1) / 2];
  else
    thickness = dist;
  return (thickness);
}

double dist_double_bonds(const vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, double avg)
{
  vector<double> a;
  int n = 0;
  double max_dist_double_bond = 0;

  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      {
        double l1 = bond_length(bond, i, atom);
        bond[i].conjoined = false;
        for (int j = i + 1; j < n_bond; j++)
          if ((bond[j].exists) && (fabs(angle_between_bonds(bond, i, j, atom)) > D_T_TOLERANCE))
            {
              double l2 = bond_length(bond, j, atom);
              double dbb = distance_between_bonds(bond, i, j, atom);

              if (dbb < avg / 2 && l1 > avg / 3 && l2 > avg / 3 && bonds_within_each_other(bond, i, j, atom))
                {
                  a.push_back(dbb);
                  n++;
                }
            }
      }
  std::sort(a.begin(), a.end());
  //for (int i = 0; i < n; i++) cout << a[i] << endl;
  //cout << "-----------------" << endl;
  if (n > 0)
    max_dist_double_bond = a[3 * (n - 1) / 4];

  if (max_dist_double_bond < 1)
    max_dist_double_bond = avg / 3;
  else
    {
      max_dist_double_bond += 2;
      for (int i = 0; i < n; i++)
        if (a[i] - max_dist_double_bond < 1 && a[i] > max_dist_double_bond)
          max_dist_double_bond = a[i];
    }
  max_dist_double_bond += 0.001;
  return (max_dist_double_bond);
}

int double_triple_bonds(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, double avg, int &n_atom,
                        double max_dist_double_bond)
{
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      {
        double l1 = bond_length(bond, i, atom);
        for (int j = i + 1; j < n_bond; j++)
          if ((bond[j].exists) && (fabs(angle_between_bonds(bond, i, j, atom)) > D_T_TOLERANCE))
            {
              double l2 = bond_length(bond, j, atom);
              double dij = distance_between_bonds(bond, i, j, atom);
              if (dij <= max_dist_double_bond && bonds_within_each_other(bond, i, j, atom))
                {
                  // start triple bond search
                  for (int k = j + 1; k < n_bond; k++)
                    if ((bond[k].exists) && (fabs(angle_between_bonds(bond, k, j, atom)) > D_T_TOLERANCE))
                      {
                        double l3 = bond_length(bond, k, atom);
                        double djk = distance_between_bonds(bond, k, j, atom);
                        double dik = distance_between_bonds(bond, k, i, atom);
                        if (djk <= max_dist_double_bond && bonds_within_each_other(bond, k, j, atom))
                          {
                            if (dik > dij)
                              {
                                bond[k].exists = false;
                                if ((l3 > l2 / 2) || (l2 > avg && l2 > 1.5 * l3 && l3 > 0.5 * avg))
                                  {
                                    bond[j].type += bond[k].type;
                                    if (bond[j].curve == bond[k].curve)
                                      bond[j].conjoined = true;
                                  }
                                if (bond[k].arom)
                                  bond[j].arom = true;
                              }
                            else
                              {
                                bond[j].exists = false;
                                if ((l2 > l3 / 2) || (l3 > avg && l3 > 1.5 * l2 && l2 > 0.5 * avg))
                                  {
                                    bond[k].type += bond[j].type;
                                    if (bond[j].curve == bond[k].curve)
                                      bond[k].conjoined = true;
                                  }
                                if (bond[j].arom)
                                  bond[k].arom = true;
                                break;
                              }
                          }
                      }

                  if (!bond[j].exists)
                    continue;
                  // end triple bond search

                  int ii = i, jj = j;
                  double l11 = l1, l22 = l2;
                  bool extended_triple = false;
                  if (l1 > avg && l1 > 1.5 * l2 && l2 > 0.5 * avg)
                    extended_triple = true;
                  else if (l2 > avg && l2 > 1.5 * l1 && l1 > 0.5 * avg)
                    {
                      ii = j;
                      jj = i;
                      l11 = l2;
                      l22 = l1;
                      extended_triple = true;
                    }
                  if (extended_triple)
                    {
                      double aa = fabs(distance_from_bond_x_a(atom[bond[ii].a].x, atom[bond[ii].a].y,
                                                              atom[bond[ii].b].x, atom[bond[ii].b].y, atom[bond[jj].a].x, atom[bond[jj].a].y));
                      double ab = fabs(distance_from_bond_x_a(atom[bond[ii].a].x, atom[bond[ii].a].y,
                                                              atom[bond[ii].b].x, atom[bond[ii].b].y, atom[bond[jj].b].x, atom[bond[jj].b].y));
                      double ba = fabs(distance_from_bond_x_b(atom[bond[ii].a].x, atom[bond[ii].a].y,
                                                              atom[bond[ii].b].x, atom[bond[ii].b].y, atom[bond[jj].a].x, atom[bond[jj].a].y));
                      double bb = fabs(distance_from_bond_x_b(atom[bond[ii].a].x, atom[bond[ii].a].y,
                                                              atom[bond[ii].b].x, atom[bond[ii].b].y, atom[bond[jj].b].x, atom[bond[jj].b].y));
                      double da = min(aa, ab);
                      double db = min(ba, bb);
                      if (da > 0.5 * l22)
                        {
                          double x = atom[bond[ii].a].x + (atom[bond[ii].b].x - atom[bond[ii].a].x) * da / l11;
                          double y = atom[bond[ii].a].y + (atom[bond[ii].b].y - atom[bond[ii].a].y) * da / l11;
			  int ba = bond[ii].a;
			  atom_t at1(x,y,bond[ii].curve);
			  at1.exists = true;
			  if (n_atom < MAX_ATOMS)
			    {
			      bond[ii].a = n_atom;
			      atom.push_back(at1);
			      n_atom++;
			    }     
			  bond_t b1(ba,bond[ii].a,bond[ii].curve);
			  if (n_bond < MAX_ATOMS)
			    {
			      bond.push_back(b1);
			      n_bond++;
			    }
                        }
                      if (db > 0.5 * l22)
                        {
                          double x = atom[bond[ii].b].x + (atom[bond[ii].a].x - atom[bond[ii].b].x) * db / l11;
                          double y = atom[bond[ii].b].y + (atom[bond[ii].a].y - atom[bond[ii].b].y) * db / l11;
			 
                          int ba = bond[ii].b;
			  atom_t at1(x,y,bond[ii].curve);
			  at1.exists = true;
			  if (n_atom < MAX_ATOMS)
			    {
			      bond[ii].b = n_atom;
			      atom.push_back(at1);
			      n_atom++;
			    }     
			  bond_t b1(ba,bond[ii].b,bond[ii].curve);
			  if (n_bond < MAX_ATOMS)
			    {
			      bond.push_back(b1);
			      n_bond++;
			    }
                        }
                      bond[jj].exists = false;
                      bond[ii].type += bond[jj].type;
                      if (bond[jj].arom)
                        bond[ii].arom = true;
                      if (bond[jj].curve == bond[ii].curve)
                        bond[ii].conjoined = true;
                      if (i == jj)
                        break;
                    }
                  else
                    {
                      if (l1 > l2)
                        {
                          bond[j].exists = false;
                          if (l2 > l1 / 2)
                            {
                              bond[i].type += bond[j].type;
                              if (bond[j].curve == bond[i].curve)
                                bond[i].conjoined = true;
                            }
                          if (bond[j].arom)
                            bond[i].arom = true;
                        }
                      else
                        {
                          bond[i].exists = false;
                          if (l1 > l2 / 2)
                            {
                              bond[j].type += bond[i].type;
                              if (bond[j].curve == bond[i].curve)
                                bond[j].conjoined = true;
                            }
                          if (bond[i].arom)
                            bond[j].arom = true;
                          break;
                        }
                    }
                }
            }
      }
  return (n_bond);
}

void extend_terminal_bond_to_label(vector<atom_t> &atom, const vector<letters_t> &letters, int n_letters, const vector<
                                   bond_t> &bond, int n_bond, const vector<label_t> &label, int n_label, double avg, double maxh,
                                   double max_dist_double_bond)
{
  for (int j = 0; j < n_bond; j++)
    if (bond[j].exists)
      {
        bool not_corner_a = terminal_bond(bond[j].a, j, bond, n_bond);
        bool not_corner_b = terminal_bond(bond[j].b, j, bond, n_bond);
        if (atom[bond[j].a].label != " ")
          not_corner_a = false;
        if (atom[bond[j].b].label != " ")
          not_corner_b = false;
        double xa = atom[bond[j].a].x;
        double ya = atom[bond[j].a].y;
        double xb = atom[bond[j].b].x;
        double yb = atom[bond[j].b].y;
        double bl = bond_length(bond, j, atom);
        double minb = FLT_MAX;
        bool found1 = false, found2 = false;
        int l1 = -1, l2 = -1;
        if (not_corner_a)
          {
            for (int i = 0; i < n_label; i++)
              if ((label[i].a)[0] != '+' && (label[i].a)[0] != '-')
                {
                  double d1 = distance_from_bond_x_a(xa, ya, xb, yb, label[i].x1, label[i].y1);
                  double d2 = distance_from_bond_x_a(xa, ya, xb, yb, label[i].x2, label[i].y2);
                  double h1 = fabs(distance_from_bond_y(xa, ya, xb, yb, label[i].x1, label[i].y1));
                  double h2 = fabs(distance_from_bond_y(xa, ya, xb, yb, label[i].x2, label[i].y2));
                  double y_dist = maxh + label[i].r1 / 2;
                  if (bond[j].type > 1)
                    y_dist += max_dist_double_bond;
                  double nb = fabs(d1) - label[i].r1 / 2;
                  if (nb <= avg && h1 <= y_dist && nb < minb && d1 < bl / 2)
                    {
                      found1 = true;
                      l1 = i;
                      minb = nb;
                    }
                  y_dist = maxh + label[i].r2 / 2;
                  if (bond[j].type > 1)
                    y_dist += max_dist_double_bond;
                  nb = fabs(d2) - label[i].r2 / 2;
                  if (nb <= avg && h2 <= y_dist && nb < minb && d2 < bl / 2)
                    {
                      found1 = true;
                      l1 = i;
                      minb = nb;
                    }
                }
            for (int i = 0; i < n_letters; i++)
              if (letters[i].free && letters[i].a != '+' && letters[i].a != '-')
                {
                  double d = distance_from_bond_x_a(xa, ya, xb, yb, letters[i].x, letters[i].y);
                  double y_dist = maxh + letters[i].r / 2;
                  if (bond[j].type > 1)
                    y_dist += max_dist_double_bond;
                  double h = fabs(distance_from_bond_y(xa, ya, xb, yb, letters[i].x, letters[i].y));
                  double nb = fabs(d) - letters[i].r / 2;
                  if (nb <= avg && h <= y_dist && nb < minb && d < bl / 2)
                    {
                      found2 = true;
                      l2 = i;
                      minb = nb;
                    }
                }
            if (found2)
              {
		//cout<< letters[l2].a<<" "<<letters[l2].x<<" "<<letters[l2].y<<" "<< atom[bond[j].a].x<<" "<< atom[bond[j].a].y<<" "<<j<<endl;
		//cout<<"=== "<<fabs(distance_from_bond_x_a(xa, ya, xb, yb, letters[l2].x, letters[l2].y))<<" "<<avg<<" "<< fabs(distance_from_bond_y(xa, ya, xb, yb, letters[l2].x, letters[l2].y))<<" "<<maxh<<" "<<letters[l2].r<<endl;
                atom[bond[j].a].label = toupper(letters[l2].a);
                atom[bond[j].a].x = letters[l2].x;
                atom[bond[j].a].y = letters[l2].y;
		atom[bond[j].a].min_x = letters[l2].min_x;
		atom[bond[j].a].min_y = letters[l2].min_y;
		atom[bond[j].a].max_x = letters[l2].max_x;
		atom[bond[j].a].max_y = letters[l2].max_y;
              }
            else if (found1)
              {
		//cout<<label[l1].a<<" "<<label[l1].x1<<" "<<label[l1].y1<<" "<<atom[bond[j].a].x<<" "<<atom[bond[j].a].y<<" "<<j<<endl;
                atom[bond[j].a].label = label[l1].a;
                atom[bond[j].a].x = (label[l1].x1 + label[l1].x2) / 2;
                atom[bond[j].a].y = (label[l1].y1 + label[l1].y2) / 2;
		atom[bond[j].a].min_x = label[l1].min_x;
		atom[bond[j].a].min_y = label[l1].min_y;
		atom[bond[j].a].max_x = label[l1].max_x;
		atom[bond[j].a].max_y = label[l1].max_y;
              }
          }
        if (not_corner_b)
          {
            found1 = false, found2 = false;
            minb = FLT_MAX;
            for (int i = 0; i < n_label; i++)
              if ((label[i].a)[0] != '+' && (label[i].a)[0] != '-' && i != l1)
                {
                  double d1 = distance_from_bond_x_b(xa, ya, xb, yb, label[i].x1, label[i].y1);
                  double d2 = distance_from_bond_x_b(xa, ya, xb, yb, label[i].x2, label[i].y2);
                  double h1 = fabs(distance_from_bond_y(xa, ya, xb, yb, label[i].x1, label[i].y1));
                  double h2 = fabs(distance_from_bond_y(xa, ya, xb, yb, label[i].x2, label[i].y2));
                  double y_dist = maxh + label[i].r1 / 2;
                  if (bond[j].type > 1)
                    y_dist += max_dist_double_bond;
                  double nb = fabs(d1) - label[i].r1 / 2; // end "b" and 1st side
                  if (nb <= avg && h1 <= y_dist && nb < minb && d1 > -bl / 2)
                    {
                      found1 = true;
                      l1 = i;
                      minb = nb;
                    }
                  y_dist = maxh + label[i].r2 / 2;
                  if (bond[j].type > 1)
                    y_dist += max_dist_double_bond;
                  nb = fabs(d2) - label[i].r2 / 2; // end "b" and 2nd side
                  if (nb <= avg && h2 <= y_dist && nb < minb && d2 > -bl / 2)
                    {
                      found1 = true;
                      l1 = i;
                      minb = nb;
                    }
                }
            for (int i = 0; i < n_letters; i++)
              if (letters[i].free && letters[i].a != '+' && letters[i].a != '-' && i != l2)
                {
                  double d = distance_from_bond_x_b(xa, ya, xb, yb, letters[i].x, letters[i].y);
                  double nb = fabs(d) - letters[i].r / 2; // distance between end "b" and letter
                  double y_dist = maxh + letters[i].r / 2;
                  if (bond[j].type > 1)
                    y_dist += max_dist_double_bond;
                  double h = fabs(distance_from_bond_y(xa, ya, xb, yb, letters[i].x, letters[i].y));
                  if (nb <= avg && h <= y_dist && nb < minb && d > -bl / 2)
                    {
                      found2 = true;
                      l2 = i;
                      minb = nb;
                    }
                }

            if (found2)
              {
		//cout<<letters[l2].a<<" "<<letters[l2].x<<" "<<letters[l2].y<<" "<< atom[bond[j].b].x<<" "<< atom[bond[j].b].y<<" "<<j<<endl;
		//cout<<"=== "<<fabs(distance_from_bond_x_a(xa, ya, xb, yb, letters[l2].x, letters[l2].y))<<" "<<avg<<" "<< fabs(distance_from_bond_y(xa, ya, xb, yb, letters[l2].x, letters[l2].y))<<" "<<maxh<<" "<<letters[l2].r<<endl;
                atom[bond[j].b].label = toupper(letters[l2].a);
                atom[bond[j].b].x = letters[l2].x;
                atom[bond[j].b].y = letters[l2].y;
		atom[bond[j].b].min_x = letters[l2].min_x;
		atom[bond[j].b].min_y = letters[l2].min_y;
		atom[bond[j].b].max_x = letters[l2].max_x;
		atom[bond[j].b].max_y = letters[l2].max_y;
              }
            else if (found1)
              {
		//cout<<label[l1].a<<" "<<label[l1].x1<<" "<<label[l1].y1<<" "<<atom[bond[j].b].x<<" "<<atom[bond[j].b].y<<" "<<j<<endl;
                atom[bond[j].b].label = label[l1].a;
                atom[bond[j].b].x = (label[l1].x1 + label[l1].x2) / 2;
                atom[bond[j].b].y = (label[l1].y1 + label[l1].y2) / 2;
		atom[bond[j].b].min_x = label[l1].min_x;
		atom[bond[j].b].min_y = label[l1].min_y;
		atom[bond[j].b].max_x = label[l1].max_x;
		atom[bond[j].b].max_y = label[l1].max_y;
              }
          }
      }
}

void extend_terminal_bond_to_bonds(vector<atom_t> &atom, vector<bond_t> &bond, int n_bond, double avg, double maxh,
                                   double max_dist_double_bond)
{
  bool found_intersection = true;

  while (found_intersection)
    {
      found_intersection = false;
      for (int j = 0; j < n_bond; j++)
        if (bond[j].exists)
          {
            bool not_corner_a = terminal_bond(bond[j].a, j, bond, n_bond);
            bool not_corner_b = terminal_bond(bond[j].b, j, bond, n_bond);
            double xa = atom[bond[j].a].x;
            double ya = atom[bond[j].a].y;
            double xb = atom[bond[j].b].x;
            double yb = atom[bond[j].b].y;
            double bl = bond_length(bond, j, atom);
            double minb = FLT_MAX;
            bool found = false;
            int l = -1;
            for (int i = 0; i < n_bond; i++)
              if (bond[i].exists && i != j)
                if (not_corner_a)
                  {
                    double h1 = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y,
                                                          atom[bond[i].b].x, atom[bond[i].b].y, xa, ya));

                    double y_dist = maxh;
                    double y_dist1 = maxh;

                    if (bond[j].type > 1 && !bond[j].conjoined)
                      y_dist += max_dist_double_bond;
                    if (bond[i].type > 1 && !bond[i].conjoined)
                      y_dist1 += max_dist_double_bond;

                    int ai = bond[i].a;
                    if (ai != bond[j].a && ai != bond[j].b)
                      {
                        double d = distance_from_bond_x_a(xa, ya, xb, yb, atom[ai].x, atom[ai].y);
                        double h = fabs(distance_from_bond_y(xa, ya, xb, yb, atom[ai].x, atom[ai].y));
                        if (fabs(d) <= avg / 2 && h <= y_dist && fabs(d) < minb && d < bl / 2 && h1 < y_dist1)
                          {
                            found = true;
                            l = ai;
                            minb = fabs(d);
                          }
                      }
                    int bi = bond[i].b;
                    if (bi != bond[j].a && bi != bond[j].b)
                      {
                        double d = distance_from_bond_x_a(xa, ya, xb, yb, atom[bi].x, atom[bi].y);
                        double h = fabs(distance_from_bond_y(xa, ya, xb, yb, atom[bi].x, atom[bi].y));
                        if (fabs(d) <= avg / 2 && h <= y_dist && fabs(d) < minb && d < bl / 2 && h1 < y_dist1)
                          {
                            found = true;
                            l = bi;
                            minb = fabs(d);
                          }
                      }
                  }
            if (found)
              {
		if (atom[l].label.length() < 2 &&  atom[bond[j].a].label.length() < 2)
		  {
		    atom[l].x = (atom[bond[j].a].x + atom[l].x) / 2;
		    atom[l].y = (atom[bond[j].a].y + atom[l].y) / 2;
		    atom[l].min_x = min(atom[bond[j].a].min_x ,atom[l].min_x);
		    atom[l].min_y = min(atom[bond[j].a].min_y ,atom[l].min_y);
		    atom[l].max_x = max(atom[bond[j].a].max_x ,atom[l].max_x);
		    atom[l].max_y = max(atom[bond[j].a].max_y ,atom[l].max_y);
		    bond[j].a = l;
		    found_intersection = true;
		  }
		//cout<<atom[bond[j].a].label<<" "<<atom[l].label<<endl;
              }

            found = false;
            minb = FLT_MAX;
            l = -1;
            for (int i = 0; i < n_bond; i++)
              if (bond[i].exists && i != j)
                if (not_corner_b)
                  {
                    double h1 = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y,
                                                          atom[bond[i].b].x, atom[bond[i].b].y, xb, yb));
                    double y_dist = maxh;
                    double y_dist1 = maxh;
                    if (bond[j].type > 1 && !bond[j].conjoined)
                      y_dist += max_dist_double_bond;
                    if (bond[i].type > 1 && !bond[i].conjoined)
                      y_dist1 += max_dist_double_bond;

                    int ai = bond[i].a;
                    if (ai != bond[j].a && ai != bond[j].b)
                      {
                        double d = distance_from_bond_x_b(xa, ya, xb, yb, atom[ai].x, atom[ai].y);
                        double h = fabs(distance_from_bond_y(xa, ya, xb, yb, atom[ai].x, atom[ai].y));
                        if (fabs(d) <= avg / 2 && h <= y_dist && fabs(d) < minb && d > -bl / 2 && h1 < y_dist1)
                          {
                            found = true;
                            l = ai;
                            minb = fabs(d);
                          }
                      }
                    int bi = bond[i].b;
                    if (bi != bond[j].a && bi != bond[j].b)
                      {
                        double d = distance_from_bond_x_b(xa, ya, xb, yb, atom[bi].x, atom[bi].y);
                        double h = fabs(distance_from_bond_y(xa, ya, xb, yb, atom[bi].x, atom[bi].y));
                        if (fabs(d) <= avg / 2 && h <= y_dist && fabs(d) < minb && d > -bl / 2 && h1 < y_dist1)
                          {
                            found = true;
                            l = bi;
                            minb = fabs(d);
                          }
                      }
                  }

            if (found)
              {
	
		if (atom[l].label.length() < 2 &&  atom[bond[j].b].label.length() < 2)
		  {
		    atom[l].x = (atom[bond[j].b].x + atom[l].x) / 2;
		    atom[l].y = (atom[bond[j].b].y + atom[l].y) / 2;
		    atom[l].min_x = min(atom[bond[j].b].min_x ,atom[l].min_x);
		    atom[l].min_y = min(atom[bond[j].b].min_y ,atom[l].min_y);
		    atom[l].max_x = max(atom[bond[j].b].max_x ,atom[l].max_x);
		    atom[l].max_y = max(atom[bond[j].b].max_y ,atom[l].max_y);
		    bond[j].b = l;
		    found_intersection = true;
		  }
		//cout<<atom[bond[j].b].label<<" "<<atom[l].label<<endl;
              }
          }
    }
}

void assign_charge(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond, const map<string, string> &fix,
                   const map<string, string> &superatom, bool debug)
{
  for (int j = 0; j < n_bond; j++)
    if (bond[j].exists && (!atom[bond[j].a].exists || !atom[bond[j].b].exists))
      bond[j].exists = false;

  for (int i = 0; i < n_atom; i++)
    if (atom[i].exists)
      {
        int n = 0;
        int m = 0;
        for (int j = 0; j < n_bond; j++)
          if (bond[j].exists && (bond[j].a == i || bond[j].b == i))
            {
              n += bond[j].type;
              if (bond[j].type > 1)
                m++;
            }
        atom[i].charge = 0;
        bool cont = true;
        while (cont)
          {
            string::size_type pos = atom[i].label.find_first_of('-');
            if (pos != string::npos)
              {
                atom[i].label.erase(pos, 1);
                if (atom[i].label.length() > 0 && isalpha(atom[i].label.at(0)))
                  atom[i].charge--;
              }
            else
              {
                pos = atom[i].label.find_first_of('+');
                if (pos != string::npos)
                  {
                    atom[i].label.erase(pos, 1);
                    if (atom[i].label.length() > 0 && isalpha(atom[i].label.at(0)))
                      atom[i].charge++;
                  }
                else
                  cont = false;
              }
          }
        for (int j = 0; j < n_bond; j++)
          if (bond[j].exists && bond[j].hash && (bond[j].b == i || bond[j].a == i))
            atom[i].charge = 0;

        atom[i].label = fix_atom_name(atom[i].label, n, fix, superatom, debug);
      }
}

int next_atom(int cur, int begin, int total)
{
  int n = cur + 1;
  if (n > total - 1)
    {
      n = begin;
    }
  return (n);
}

bool dir_change(int n, int last, int begin, int total, const vector<atom_t> &atom)
{
  int m = next_atom(n, begin, total);
  while (distance(atom[m].x, atom[m].y, atom[n].x, atom[n].y) < V_DISPLACEMENT && m != n)
    m = next_atom(m, begin, total);
  if (m == n)
    return (false);
  double s = fabs(distance_from_bond_y(atom[n].x, atom[n].y, atom[last].x, atom[last].y, atom[m].x, atom[m].y));
  if (s > DIR_CHANGE)
    return (true);
  return (false);
}

bool smaller_distance(int n, int last, int begin, int total, const vector<atom_t> &atom)
{
  int m = next_atom(n, begin, total);
  double d1 = distance(atom[n].x, atom[n].y, atom[last].x, atom[last].y);
  double d2 = distance(atom[m].x, atom[m].y, atom[last].x, atom[last].y);
  if (d1 > d2)
    {
      return (true);
    }
  return (false);
}

void douglas_peucker(vector<atom_t> &atom,  int b_atom, int n_atom)
{
  double max_s = 0;
  double max_i = b_atom;
  for (int i=b_atom+1; i<n_atom; i++)
    {
      double s = fabs(distance_from_bond_y(atom[b_atom].x, atom[b_atom].y, atom[n_atom-1].x, atom[n_atom-1].y, atom[i].x, atom[i].y));
      if (s>max_s)
	{
	  max_i = i;
	  max_s = s;
	}
    }
  if (max_s>DIR_CHANGE)
    {
      vector<atom_t> atom1(atom);
      vector<atom_t> atom2(atom);
      douglas_peucker(atom1,  b_atom, max_i);
      douglas_peucker(atom2,  max_i, n_atom);
      for (int i=b_atom; i<max_i; i++)
	atom[i].exists = atom1[i].exists;
      for (int i=max_i; i<n_atom; i++)
	atom[i].exists = atom2[i].exists;
    }
  else
    {
      atom[b_atom].exists = true;
      atom[n_atom-1].exists = true;
    }
}


int find_bonds(vector<atom_t> &atom, vector<bond_t> &bond, int b_atom, int n_atom, int n_bond,
               const potrace_path_t * const p)
{

  //douglas_peucker(atom,  b_atom, n_atom);
  
  int i = b_atom + 1;
  int last = b_atom;

  while (i < n_atom)
    {
      if (atom[i].corner)
        {
          atom[i].exists = true;
          last = i;
          i++;
        }
      else if (dir_change(i, last, b_atom, n_atom, atom))
        {
          atom[i].exists = true;
          last = i;
          i++;
        }
      else if (smaller_distance(i, last, b_atom, n_atom, atom))
        {
          atom[i].exists = true;
          last = i;
          i++;
        }
      else
        {
          i++;
        }
    }
  
  for (int i = b_atom; i < n_atom; i++)
    if (atom[i].exists)
      {
	int j = next_atom(i, b_atom, n_atom);
        while (!atom[j].exists)
          {
            j = next_atom(j, b_atom, n_atom);
          }
        bond_t bn(i,j,p);
        if (n_bond < MAX_ATOMS)
	  {
	    bond.push_back(bn);
	    n_bond++;
	  }
      }
  return (n_bond);
}



int find_atoms(const potrace_path_t *p, vector<atom_t> &atom, vector<bond_t> &bond, int *n_bond,int width, int height)
{
  int *tag, n_atom = 0;
  potrace_dpoint_t (*c)[3];
  long n;

  while (p != NULL)
    {
      n = p->curve.n;
      tag = p->curve.tag;
      c = p->curve.c;
      int b_atom = n_atom;
      double x = c[n - 1][2].x;
      double y = c[n - 1][2].y;
      if (x<0) x=0;
      if (x>width) x=width;
      if (y<0) y=0;
      if (y>height) y=height;
      atom_t at(x,y,p);

      if (n_atom < MAX_ATOMS)
	{
	  atom.push_back(at);
	  n_atom++;
	}     

      for (long i = 0; i < n; i++)
        {

          switch (tag[i])
            {
            case POTRACE_CORNER:
	      {
		x = c[i][1].x;
		y = c[i][1].y;
		if (x<0) x=0;
		if (x>width) x=width;
		if (y<0) y=0;
		if (y>height) y=height;
		atom_t at1(x,y,p);
		at1.corner = true;
		if (n_atom < MAX_ATOMS)
		  {
		    atom.push_back(at1);
		    n_atom++;
		  }     
	      }
              break;
            case POTRACE_CURVETO:
	      {
		x = c[i][0].x;
		y = c[i][0].y;
		if (x<0) x=0;
		if (x>width) x=width;
		if (y<0) y=0;
		if (y>height) y=height;
		atom_t at2(x,y,p);

		if (n_atom < MAX_ATOMS)
		  {
		    atom.push_back(at2);
		    n_atom++;
		  }     
		
		x = c[i][1].x;
		y = c[i][1].y;
		if (x<0) x=0;
		if (x>width) x=width;
		if (y<0) y=0;
		if (y>height) y=height;
		atom_t at3(x,y,p);
		
		if (n_atom < MAX_ATOMS)
		  {
		    atom.push_back(at3);
		    n_atom++;
		  }     
	      }
              break;
            }
          if (i != n - 1)
            {
	      x = c[i][2].x;
	      y = c[i][2].y;
	      if (x<0) x=0;
	      if (x>width) x=width;
	      if (y<0) y=0;
	      if (y>height) y=height;
	      atom_t at4(x,y,p);

	      if (n_atom < MAX_ATOMS)
		{
		  atom.push_back(at4);
		  n_atom++;
		}     
            }
        }
      *n_bond = find_bonds(atom, bond, b_atom, n_atom, *n_bond, p);
      p = p->next;
    }
  return (n_atom);
}

bool comp_dashes_x(const dash_t &a, const dash_t &b)
{
  return (a.x < b.x);
}

bool comp_dashes_y(const dash_t &a, const dash_t &b)
{
  return (a.y < b.y);
}

void extend_dashed_bond(int a, int b, int n, vector<atom_t> &atom)
{
  double x0 = atom[a].x;
  double y0 = atom[a].y;
  double x1 = atom[b].x;
  double y1 = atom[b].y;
  double l = distance(x0, y0, x1, y1);
  double kx = (x1 - x0) / l;
  double ky = (y1 - y0) / l;
  double x = kx * (-1. * l / (n - 1)) + x0;
  double y = ky * (-1. * l / (n - 1)) + y0;
  atom[a].x = x;
  atom[a].y = y;
  atom[b].x = kx * l / (n - 1) + x1;
  atom[b].y = ky * l / (n - 1) + y1;
}

int count_area(vector<vector<int> > &box, double &x0, double &y0)
{
  int a = 0;
  int w = box.size();
  int h = box[0].size();
  int x = int(x0);
  int y = int(y0);
  int xm = 0, ym = 0;

  if (box[x][y] == 1)
    {
      box[x][y] = 2;
      list<int> cx;
      list<int> cy;
      cx.push_back(x);
      cy.push_back(y);
      while (!cx.empty())
        {
          x = cx.front();
          y = cy.front();
          cx.pop_front();
          cy.pop_front();
          box[x][y] = 0;
          a++;
          xm += x;
          ym += y;
          for (int i = x - 1; i < x + 2; i++)
            for (int j = y - 1; j < y + 2; j++)
              if (i < w && j < h && i >= 0 && j >= 0 && box[i][j] == 1)
                {
                  cx.push_back(i);
                  cy.push_back(j);
                  box[i][j] = 2;
                }
        }
    }
  else
    return (0);

  x0 = 1. * xm / a;
  y0 = 1. * ym / a;

  return (a);
}

bool curve_in_letter(const potrace_path_t *p,  vector<letters_t> &letters)
{
  bool r = false;
  for (int i=0; i<letters.size(); i++)
    if (letters[i].curve == p)
      r = true;
  return r;
}

int find_dashed_bonds(const potrace_path_t *p, vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int *n_bond,
                      int max, double avg, const Image &img, const ColorGray &bg, double THRESHOLD, bool thick, double dist, vector<letters_t> &letters)
{
  potrace_dpoint_t (*c)[3];
  vector<dash_t> dot;
  vector<vector<int> > box(img.columns());
  int width = img.columns();
  int height = img.rows();

  for (unsigned int i = 0; i < img.columns(); i++)
    for (unsigned int j = 0; j < img.rows(); j++)
      box[i].push_back(get_pixel(img, bg, i, j, THRESHOLD));

  while (p != NULL)
    {
      if (p->sign == int('+') && p->area < max && p->childlist == NULL && !curve_in_letter(p,letters) && dot.size() < 100)
        {
	  int n = p->curve.n;
          c = p->curve.c;
          int *tag = p->curve.tag;
          double cx =  c[n - 1][2].x;
          double cy =  c[n - 1][2].y;
          if (cx<0) cx=0;
          if (cx>width) cx=width;
          if (cy<0) cy=0;
          if (cy>height) cy=height;
	  dash_t d;
          d.x = cx;
          d.y = cy;
          double l = cx;
          double r = cx;
          double t = cy;
          double b = cy;
          d.curve = p;
          d.free = true;

          int tot = 1;
          for (long i = 0; i < n; i++)
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
                  d.x += cx;
                  d.y += cy;
                  if (cx < l)
                    l = cx;
                  if (cx > r)
                    r = cx;
                  if (cy < t)
                    t = cy;
                  if (cx > b)
                    b = cy;
                  tot++;
                  break;
                case POTRACE_CURVETO:
                  cx = c[i][0].x;
                  cy = c[i][0].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  d.x += cx;
                  d.y += cy;
                  if (cx < l)
                    l = cx;
                  if (cx > r)
                    r = cx;
                  if (cy < t)
                    t = cy;
                  if (cx > b)
                    b = cy;
                  cx = c[i][1].x;
                  cy = c[i][1].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  d.x += cx;
                  d.y += cy;
                  if (cx < l)
                    l = cx;
                  if (cx > r)
                    r = cx;
                  if (cy < t)
                    t = cy;
                  if (cx > b)
                    b = cy;
                  tot += 2;
                  break;
                }
              if (i != n - 1)
                {
                  cx = c[i][2].x;
                  cy = c[i][2].y;
                  if (cx<0) cx=0;
                  if (cx>width) cx=width;
                  if (cy<0) cy=0;
                  if (cy>height) cy=height;
                  d.x += cx;
                  d.y += cy;
                  if (cx < l)
                    l = cx;
                  if (cx > r)
                    r = cx;
                  if (cy < t)
                    t = cy;
                  if (cx > b)
                    b = cy;
                  tot++;
                }
            }
          d.x /= tot;
          d.y /= tot;
          if (thick)
            d.area = count_area(box, d.x, d.y);
          else
            d.area = p->area;
          if (distance(l, t, r, b) < avg / 3)
	    dot.push_back(d);
        }
      p = p->next;
    }
  for (int i = 0; i < dot.size(); i++)
    if (dot[i].free)
      {
        vector<dash_t> dash;
        dash.push_back(dot[i]);
        dot[i].free = false;
        double l = dot[i].x;
        double r = dot[i].x;
        double t = dot[i].y;
        double b = dot[i].y;
        double mx = l;
        double my = t;
        double dist_next = FLT_MAX;
        int next_dot = i;
        for (int j = i + 1; j < dot.size(); j++)
          if (dot[j].free && distance(dash[0].x, dash[0].y, dot[j].x, dot[j].y) <= dist && distance(dash[0].x,
              dash[0].y, dot[j].x, dot[j].y) < dist_next)
            {
              dash.push_back(dot[j]);
              dist_next = distance(dash[0].x, dash[0].y, dot[j].x, dot[j].y);
              next_dot = j;
            }

        if (next_dot != i)
          {
            dot[next_dot].free = false;
            if (dash[1].x < l)
              l = dash[1].x;
            if (dash[1].x > r)
              r = dash[1].x;
            if (dash[1].y < t)
              t = dash[1].y;
            if (dash[1].y > b)
              b = dash[1].y;
            mx = (mx + dash[1].x) / 2;
            my = (my + dash[1].y) / 2;
          }
        bool found = true;
        while (dash.size() > 1 && found)
          {
            dist_next = FLT_MAX;
            found = false;
            int minj = next_dot;
	    dash_t d;
            for (int j = next_dot + 1; j < dot.size(); j++)
              if (dot[j].free && distance(mx, my, dot[j].x, dot[j].y) <= dist && distance(mx, my, dot[j].x,
                  dot[j].y) < dist_next
                  && fabs(distance_from_bond_y(dash[0].x, dash[0].y, dash.back().x, dash.back().y, dot[j].x,
                                               dot[j].y)) < V_DISPLACEMENT)
                {
                  d = dot[j];
                  dist_next = distance(mx, my, dot[j].x, dot[j].y);
                  found = true;
                  minj = j;
                }
            if (found)
              {
                dot[minj].free = false;
                if (d.x < l)
                  l = d.x;
                if (d.x > r)
                  r = d.x;
                if (d.y < t)
                  t = d.y;
                if (d.y > b)
                  b = d.y;
                mx = (mx + d.x) / 2;
                my = (my + d.y) / 2;
                dash.push_back(d);
              }
          }

        if (dash.size() > 2)
          {
            if ((r - l) > (b - t))
              {
                sort(dash.begin(),dash.end(),comp_dashes_x);
              }
            else
              {
                sort(dash.begin(),dash.end(), comp_dashes_y);
              }
            bool one_line = true;
            double dx = dash.back().x - dash[0].x;
            double dy = dash.back().y - dash[0].y;
            double k = 0;
            if (fabs(dx) > fabs(dy))
              k = dy / dx;
            else
              k = dx / dy;
            for (int j = 1; j < dash.size() - 1; j++)
              {
                double nx = dash[j].x - dash[0].x;
                double ny = dash[j].y - dash[0].y;
                double diff = 0;
                if (fabs(dx) > fabs(dy))
                  diff = k * nx - ny;
                else
                  diff = k * ny - nx;
                if (fabs(diff) > V_DISPLACEMENT || dash[j].curve == dash[0].curve)
                  one_line = false;
              }
            if (one_line)
              {
                for (int j = 0; j < dash.size(); j++)
                  delete_curve(atom, bond, n_atom, *n_bond, dash[j].curve);
		atom_t at1(dash[0].x,dash[0].y,dash[0].curve);
		at1.exists = true;
		if (n_atom < MAX_ATOMS)
		  {
		    atom.push_back(at1);
		    n_atom++;
		  }     
		atom_t at2(dash.back().x,dash.back().y,dash.back().curve);
		at2.exists = true;
		if (n_atom < MAX_ATOMS)
		  {
		    atom.push_back(at2);
		    n_atom++;
		  }     
		bond_t b1(n_atom-2,n_atom-1,dash[0].curve);
		b1.hash = true;                        
		if (*n_bond < MAX_ATOMS)
		  {
		    bond.push_back(b1);
		    if (dash[0].area > dash.back().area)
		      bond_end_swap(bond, *n_bond);
		    extend_dashed_bond(bond[*n_bond].a, bond[*n_bond].b, dash.size(), atom);
		    (*n_bond)++;
		  }
              }
          }
      }

  return (n_atom);
}


int find_small_bonds(const potrace_path_t *p, vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int *n_bond,
                     double max_area, double Small, double thickness)
{
  while (p != NULL)
    {
      if ((p->sign == int('+')) && (p->area <= max_area))
        {
          vector<dash_t> dot;
          for (int i = 0; i < n_atom; i++)
            if ((atom[i].exists) && (atom[i].curve == p) && (dot.size() < 20))
              {
		dash_t d;
                d.x = atom[i].x;
                d.y = atom[i].y;
                d.curve = p;
                d.free = true;
		dot.push_back(d);
              }

          if ((dot.size() > 2))
            {
              double l = dot[0].x;
              double r = dot[0].x;
              double t = dot[0].y;
              double b = dot[0].y;
              for (int i = 1; i < dot.size(); i++)
                {
                  if (dot[i].x < l)
                    l = dot[i].x;
                  if (dot[i].x > r)
                    r = dot[i].x;
                  if (dot[i].y < t)
                    t = dot[i].y;
                  if (dot[i].y > b)
                    b = dot[i].y;
                }
              if ((r - l) > (b - t))
                {
                  sort(dot.begin(),dot.end(),comp_dashes_x);
                }
              else
                {
                  sort(dot.begin(),dot.end(), comp_dashes_y);
                }
              double d = 0;
              for (int i = 1; i < dot.size() - 1; i++)
                d = max(d, fabs(distance_from_bond_y(dot[0].x, dot[0].y, dot.back().x, dot.back().y,
                                                     dot[i].x, dot[i].y)));
              if (d < thickness || p->area < Small)
                {
                  delete_curve(atom, bond, n_atom, *n_bond, p);
		  atom_t at1(dot[0].x,dot[0].y,p);
		  at1.exists = true;
		  if (n_atom < MAX_ATOMS)
		    {
		      atom.push_back(at1);
		      n_atom++;
		    }     
		  atom_t at2(dot.back().x,dot.back().y,p);
		  at2.exists = true;
		  if (n_atom < MAX_ATOMS)
		    {
		      atom.push_back(at2);
		      n_atom++;
		    }     
		  bond_t b1(n_atom-2,n_atom-1,p);
		  b1.Small = true;
		  if (*n_bond < MAX_ATOMS)
		    {
		      bond.push_back(b1);
		      (*n_bond)++;
		    }
                }
            }
        }
      p = p->next;
    }
  return (n_atom);
}

int resolve_bridge_bonds(vector<atom_t> &atom, int n_atom, vector<bond_t> &bond, int n_bond, double thickness,
                         double avg_bond_length, const map<string, string> &superatom, bool verbose)
{
  molecule_statistics_t molecule_statistics1 = calculate_molecule_statistics(atom, bond, n_bond, avg_bond_length, superatom, verbose);

  for (int i = 0; i < n_atom; i++)
    if ((atom[i].exists) && (atom[i].label == " "))
      {
        list<int> con;
        for (int j = 0; j < n_bond; j++)
          if ((bond[j].exists) && (bond[j].a == i || bond[j].b == i))
            con.push_back(j);
        if (con.size() == 4)
          {
            int a = con.front();
            con.pop_front();
            int b = 0;
            int e = 0;
            while ((con.size() > 2) && (e++ < 3))
              {
                b = con.front();
                con.pop_front();
                double y1 = distance_from_bond_y(atom[bond[a].a].x, atom[bond[a].a].y, atom[bond[a].b].x,
                                                 atom[bond[a].b].y, atom[bond[b].a].x, atom[bond[b].a].y);
                double y2 = distance_from_bond_y(atom[bond[a].a].x, atom[bond[a].a].y, atom[bond[a].b].x,
                                                 atom[bond[a].b].y, atom[bond[b].b].x, atom[bond[b].b].y);
                if (fabs(y1) > thickness || fabs(y2) > thickness)
                  con.push_back(b);
              }
            if (con.size() == 2)
              {
                int c = con.front();
                con.pop_front();
                int d = con.front();
                con.pop_front();
                vector<int> term;
                term.push_back(a);
                term.push_back(b);
                term.push_back(c);
                term.push_back(d);
                bool terminal = false;
                for (unsigned int k = 0; k < term.size(); k++)
                  {
                    bool terminal_a = terminal_bond(bond[term[k]].a, term[k], bond, n_bond);
                    bool terminal_b = terminal_bond(bond[term[k]].b, term[k], bond, n_bond);
                    if (terminal_a || terminal_b)
                      terminal = true;
                  }
                double y1 = distance_from_bond_y(atom[bond[c].a].x, atom[bond[c].a].y, atom[bond[c].b].x,
                                                 atom[bond[c].b].y, atom[bond[d].a].x, atom[bond[d].a].y);
                double y2 = distance_from_bond_y(atom[bond[c].a].x, atom[bond[c].a].y, atom[bond[c].b].x,
                                                 atom[bond[c].b].y, atom[bond[d].b].x, atom[bond[d].b].y);
                if (bond[a].type == 1 && bond[b].type == 1 && bond[c].type == 1 && bond[d].type == 1 && fabs(y1)
                    < thickness && fabs(y2) < thickness && !terminal)
                  {
                    bond[b].exists = false;
                    bond[d].exists = false;
                    atom[i].exists = false;
                    if (bond[a].a == bond[b].a)
                      bond[a].a = bond[b].b;
                    else if (bond[a].a == bond[b].b)
                      bond[a].a = bond[b].a;
                    else if (bond[a].b == bond[b].a)
                      bond[a].b = bond[b].b;
                    else if (bond[a].b == bond[b].b)
                      bond[a].b = bond[b].a;
                    if (bond[c].a == bond[d].a)
                      bond[c].a = bond[d].b;
                    else if (bond[c].a == bond[d].b)
                      bond[c].a = bond[d].a;
                    else if (bond[c].b == bond[d].a)
                      bond[c].b = bond[d].b;
                    else if (bond[c].b == bond[d].b)
                      bond[c].b = bond[d].a;

                    molecule_statistics_t molecule_statistics2 = calculate_molecule_statistics(atom, bond, n_bond, avg_bond_length, superatom, verbose);
                    if (molecule_statistics1.fragments != molecule_statistics2.fragments ||
                        molecule_statistics1.rotors != molecule_statistics2.rotors ||
                        molecule_statistics1.rings56 - molecule_statistics2.rings56 == 2)
                      {
                        bond[b].exists = true;
                        bond[d].exists = true;
                        atom[i].exists = true;
                        if (bond[a].a == bond[b].a)
                          bond[a].a = bond[b].b;
                        else if (bond[a].a == bond[b].b)
                          bond[a].a = bond[b].a;
                        else if (bond[a].b == bond[b].a)
                          bond[a].b = bond[b].b;
                        else if (bond[a].b == bond[b].b)
                          bond[a].b = bond[b].a;
                        if (bond[c].a == bond[d].a)
                          bond[c].a = bond[d].b;
                        else if (bond[c].a == bond[d].b)
                          bond[c].a = bond[d].a;
                        else if (bond[c].b == bond[d].a)
                          bond[c].b = bond[d].b;
                        else if (bond[c].b == bond[d].b)
                          bond[c].b = bond[d].a;
                      }
                  }
              }
          }
      }
  return (molecule_statistics1.fragments);
}

void collapse_atoms(vector<atom_t> &atom, vector<bond_t> &bond, int n_atom, int n_bond, double dist)
{
  bool found = true;

  while (found)
    {
      found = false;
      for (int i = 0; i < n_atom; i++)
        if (atom[i].exists)
          for (int j = 0; j < n_atom; j++)
            if (atom[j].exists && j != i && distance(atom[i].x, atom[i].y, atom[j].x, atom[j].y) < dist)
              {
                atom[j].exists = false;
                atom[i].x = (atom[i].x + atom[j].x) / 2;
                atom[i].y = (atom[i].y + atom[j].y) / 2;

		atom[i].min_x = min(atom[i].min_x,atom[j].min_x);
		atom[i].min_y = min(atom[i].min_y,atom[j].min_y);
		atom[i].max_x = max(atom[i].max_x,atom[j].max_x);
		atom[i].max_y = max(atom[i].max_y,atom[j].max_y);

                if (atom[j].label != " " && atom[i].label == " ")
                  atom[i].label = atom[j].label;
                for (int k = 0; k < n_bond; k++)
                  if (bond[k].exists)
                    {
                      if (bond[k].a == j)
                        {
                          bond[k].a = i;
                        }
                      else if (bond[k].b == j)
                        {
                          bond[k].b = i;
                        }
                    }
                found = true;
              }
    }
}

void collapse_bonds(vector<atom_t> &atom, const vector<bond_t> &bond, int n_bond, double dist)
{
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && bond_length(bond, i, atom) < dist)
      {
        atom[bond[i].a].x = (atom[bond[i].a].x + atom[bond[i].b].x) / 2;
        atom[bond[i].a].y = (atom[bond[i].a].y + atom[bond[i].b].y) / 2;
	atom[bond[i].a].min_x = min(atom[bond[i].a].min_x ,atom[bond[i].b].min_x) ;
        atom[bond[i].a].min_y = min(atom[bond[i].a].min_y,atom[bond[i].b].min_y);
	atom[bond[i].a].max_x = max(atom[bond[i].a].max_x ,atom[bond[i].b].max_x) ;
        atom[bond[i].a].max_y = max(atom[bond[i].a].max_y,atom[bond[i].b].max_y);
        atom[bond[i].b].x = (atom[bond[i].a].x + atom[bond[i].b].x) / 2;
        atom[bond[i].b].y = (atom[bond[i].a].y + atom[bond[i].b].y) / 2;
	atom[bond[i].b].min_x = min(atom[bond[i].a].min_x ,atom[bond[i].b].min_x) ;
        atom[bond[i].b].min_y = min(atom[bond[i].a].min_y,atom[bond[i].b].min_y);
	atom[bond[i].b].max_x = max(atom[bond[i].a].max_x ,atom[bond[i].b].max_x) ;
        atom[bond[i].b].max_y = max(atom[bond[i].a].max_y,atom[bond[i].b].max_y);
      }
}


int fix_one_sided_bonds(vector<bond_t> &bond, int n_bond, const vector<atom_t> &atom, double thickness, double avg)
{
  double l;
  int n_bond_orig = n_bond;
  for (int i = 0; i < n_bond_orig; i++)
    if (bond[i].exists && bond[i].type < 3 && (l = bond_length(bond, i, atom)) > avg / 3)
      for (int j = 0; j < n_bond_orig; j++)
        if (bond[j].exists && j != i && bond[j].type < 3 && fabs(angle_between_bonds(bond, i, j, atom))
            < D_T_TOLERANCE && bond_length(bond, j, atom) > avg / 3)
          {
            double d1 = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                  atom[bond[i].b].y, atom[bond[j].a].x, atom[bond[j].a].y));
            double d2 = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                  atom[bond[i].b].y, atom[bond[j].b].x, atom[bond[j].b].y));
            if (d1 < thickness && !(bond[j].a == bond[i].b || bond[j].a == bond[i].a))
              {
                double l1 = distance_from_bond_x_a(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                   atom[bond[i].b].y, atom[bond[j].a].x, atom[bond[j].a].y);
                if (l1 > 0 && l1 < l)
                  {
                    if (bond[j].b == bond[i].b || bond[j].b == bond[i].a)
                      {
                        bond[j].exists = false;
                      }
                    else
                      {
			bond_t b1(bond[j].a,bond[i].b,bond[i].curve);
                        b1.type = bond[i].type;
			b1.hash = bond[i].hash;                        
			b1.wedge = bond[i].wedge;
			b1.arom = bond[i].arom;
                        b1.conjoined = bond[i].conjoined;
			if (n_bond < MAX_ATOMS)
			  {
			    bond.push_back(b1);
			    n_bond++;
			  }
                        bond[i].b = bond[j].a;
                        bond[i].wedge = false;
                      }
                  }
              }
            else if (d2 < thickness && !(bond[j].b == bond[i].b || bond[j].b == bond[i].a))
              {
                double l1 = distance_from_bond_x_a(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                   atom[bond[i].b].y, atom[bond[j].b].x, atom[bond[j].b].y);
                if (l1 > 0 && l1 < l)
                  {
                    if (bond[j].a == bond[i].b || bond[j].a == bond[i].a)
                      {
                        bond[j].exists = false;
                      }
                    else
                      {
			bond_t b1(bond[j].b,bond[i].b,bond[i].curve);
                        b1.type = bond[i].type;
			b1.hash = bond[i].hash;                        
			b1.wedge = bond[i].wedge;
			b1.arom = bond[i].arom;
                        b1.conjoined = bond[i].conjoined;
			if (n_bond < MAX_ATOMS)
			  {
			    bond.push_back(b1);
			    n_bond++;
			  }
                        bond[i].b = bond[j].b;
                        bond[i].wedge = false;
                      }
                  }
              }
          }

  return (n_bond);
}

int thickness_hor(const Image &image, int x1, int y1, const ColorGray &bgColor, double THRESHOLD_BOND)
{
  int i = 0, s = 0, w = 0;
  int width = image.columns();
  s = get_pixel(image, bgColor, x1, y1, THRESHOLD_BOND);

  if (s == 0 && x1 + 1 < width)
    {
      x1++;
      s = get_pixel(image, bgColor, x1, y1, THRESHOLD_BOND);
    }
  if (s == 0 && x1 - 2 >= 0)
    {
      x1 -= 2;
      s = get_pixel(image, bgColor, x1, y1, THRESHOLD_BOND);
    }
  if (s == 1)
    {
      while (x1 + i < width && s == 1)
        s = get_pixel(image, bgColor, x1 + i++, y1, THRESHOLD_BOND);
      w = i - 1;
      i = 1;
      s = 1;
      while (x1 - i >= 0 && s == 1)
        s = get_pixel(image, bgColor, x1 - i++, y1, THRESHOLD_BOND);
      w += i - 1;
    }
  return (w);
}

int thickness_ver(const Image &image, int x1, int y1, const ColorGray &bgColor, double THRESHOLD_BOND)
{
  int i = 0, s = 0, w = 0;
  int height = image.rows();
  s = get_pixel(image, bgColor, x1, y1, THRESHOLD_BOND);

  if (s == 0 && y1 + 1 < height)
    {
      y1++;
      s = get_pixel(image, bgColor, x1, y1, THRESHOLD_BOND);
    }
  if (s == 0 && y1 - 2 >= 0)
    {
      y1 -= 2;
      s = get_pixel(image, bgColor, x1, y1, THRESHOLD_BOND);
    }
  if (s == 1)
    {
      while (y1 + i < height && s == 1)
        s = get_pixel(image, bgColor, x1, y1 + i++, THRESHOLD_BOND);
      w = i - 1;
      i = 1;
      s = 1;
      while (y1 - i >= 0 && s == 1)
        s = get_pixel(image, bgColor, x1, y1 - i++, THRESHOLD_BOND);
      w += i - 1;
    }
  return (w);
}

bool comp_thickness(const pair<int,double> &a, const pair<int,double> &b)
{
  return a.second < b.second;
}

double find_wedge_bonds(const Image &image, vector<atom_t> &atom, int n_atom, vector<bond_t> &bond, int n_bond,
                        const ColorGray &bgColor, double THRESHOLD_BOND, double max_dist_double_bond, double avg, int limit, int dist)
{
  double l;
  vector<pair<int,double> > a;
  int n = 0;
  a.push_back(make_pair(-1,1.5));
  vector<int> x_reg, y_reg;
  bool triangle = false;
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && !bond[i].hash && bond[i].type == 1 && (l = bond_length(bond, i, atom))
        > max_dist_double_bond)
      {
        x_reg.clear();
        y_reg.clear();
        double avg_x = 0, avg_y = 0;
        int x1 = int((atom[bond[i].a].x + atom[bond[i].b].x) / 2);
        int y1 = int((atom[bond[i].a].y + atom[bond[i].b].y) / 2);

        int w = 0, max_c, min_c, sign = 1;
        int w3_ver = thickness_ver(image, x1, y1, bgColor, THRESHOLD_BOND);
        int w3_hor = thickness_hor(image, x1, y1, bgColor, THRESHOLD_BOND);
        if (w3_ver == 0 && w3_hor == 0)
          continue;
        if ((w3_ver < w3_hor && w3_ver > 0) || w3_hor == 0)
          {
            w = w3_ver;
            int old = w3_ver;
            max_c = int(max(atom[bond[i].a].x, atom[bond[i].b].x)) - dist;
            min_c = int(min(atom[bond[i].a].x, atom[bond[i].b].x)) + dist;
            if (atom[bond[i].b].x < atom[bond[i].a].x)
              sign = -1;
            for (int j = x1 + 1; j <= max_c; j++)
              {
                int y = int(atom[bond[i].a].y + (atom[bond[i].b].y - atom[bond[i].a].y) * (j - atom[bond[i].a].x)
                            / (atom[bond[i].b].x - atom[bond[i].a].x));
                int t = thickness_ver(image, j, y, bgColor, THRESHOLD_BOND);
                if (abs(t - old) > 2)
                  break;
                if (t < 2 * MAX_BOND_THICKNESS && t < avg / 3 && t > 0)
                  {
                    x_reg.push_back(j);
                    y_reg.push_back(t);
                    avg_x += j;
                    avg_y += t;
                    w = max(w, t);
                  }
                old = t;
              }
            old = w3_ver;
            for (int j = x1 - 1; j >= min_c; j--)
              {
                int y = int(atom[bond[i].a].y + (atom[bond[i].b].y - atom[bond[i].a].y) * (j - atom[bond[i].a].x)
                            / (atom[bond[i].b].x - atom[bond[i].a].x));
                int t = thickness_ver(image, j, y, bgColor, THRESHOLD_BOND);
                if (abs(t - old) > 2)
                  break;
                if (t < 2 * MAX_BOND_THICKNESS && t < avg / 3 && t > 0)
                  {
                    x_reg.push_back(j);
                    y_reg.push_back(t);
                    avg_x += j;
                    avg_y += t;
                    w = max(w, t);
                  }
                old = t;
              }

          }
        else
          {
            w = w3_hor;
            int old = w3_hor;
            max_c = int(max(atom[bond[i].a].y, atom[bond[i].b].y)) - dist;
            min_c = int(min(atom[bond[i].a].y, atom[bond[i].b].y)) + dist;
            if (atom[bond[i].b].y < atom[bond[i].a].y)
              sign = -1;
            for (int j = y1 + 1; j <= max_c; j++)
              {
                int x = int(atom[bond[i].a].x + (atom[bond[i].b].x - atom[bond[i].a].x) * (j - atom[bond[i].a].y)
                            / (atom[bond[i].b].y - atom[bond[i].a].y));
                int t = thickness_hor(image, x, j, bgColor, THRESHOLD_BOND);
                if (abs(t - old) > 2)
                  break;
                if (t < 2 * MAX_BOND_THICKNESS && t < avg / 3 && t > 0)
                  {
                    x_reg.push_back(j);
                    y_reg.push_back(t);
                    avg_x += j;
                    avg_y += t;
                    w = max(w, t);
                  }
                old = t;
              }
            old = w3_hor;
            for (int j = y1 - 1; j >= min_c; j--)
              {
                int x = int(atom[bond[i].a].x + (atom[bond[i].b].x - atom[bond[i].a].x) * (j - atom[bond[i].a].y)
                            / (atom[bond[i].b].y - atom[bond[i].a].y));
                int t = thickness_hor(image, x, j, bgColor, THRESHOLD_BOND);
                if (abs(t - old) > 2)
                  break;
                if (t < 2 * MAX_BOND_THICKNESS && t < avg / 3 && t > 0)
                  {
                    x_reg.push_back(j);
                    y_reg.push_back(t);
                    avg_x += j;
                    avg_y += t;
                    w = max(w, t);
                  }
                old = t;
              }
          }
        avg_x /= x_reg.size();
        avg_y /= y_reg.size();
        double numerator = 0, denominator = 0;
        for (unsigned int j = 0; j < x_reg.size(); j++)
          {
            numerator += 1. * (x_reg[j] - avg_x) * (y_reg[j] - avg_y);
            denominator += 1. * (x_reg[j] - avg_x) * (x_reg[j] - avg_x);
          }
        double beta = 0;
        if (denominator != 0)
          beta = numerator / denominator;
        //cout << fabs(beta) * (max_c - min_c) << " " << (max_c - min_c) << " " << avg << endl;
        if (fabs(beta) * (max_c - min_c) > limit)
          {
            bond[i].wedge = true;
	    triangle = true;
            if (beta * sign < 0)
              bond_end_swap(bond, i);
          }

        if (bond[i].wedge)
          {
            for (int j = 0; j < n_atom; j++)
              if (atom[j].exists && j != bond[i].b && distance(atom[bond[i].b].x, atom[bond[i].b].y, atom[j].x,
                  atom[j].y) <= w)
                {
                  atom[j].exists = false;
                  atom[bond[i].b].x = (atom[bond[i].b].x + atom[j].x) / 2;
                  atom[bond[i].b].y = (atom[bond[i].b].y + atom[j].y) / 2;
		  atom[bond[i].b].min_x = min(atom[bond[i].b].min_x, atom[j].min_x);
                  atom[bond[i].b].min_y = min(atom[bond[i].b].min_y,atom[j].min_y);
		  atom[bond[i].b].max_x = max(atom[bond[i].b].max_x, atom[j].max_x);
		  atom[bond[i].b].max_y = max(atom[bond[i].b].max_y,atom[j].max_y);
                  for (int k = 0; k < n_bond; k++)
                    if (bond[k].exists)
                      {
                        if (bond[k].a == j)
                          {
                            bond[k].a = bond[i].b;
                          }
                        else if (bond[k].b == j)
                          {
                            bond[k].b = bond[i].b;
                          }
                      }
                }
          }
        if (!bond[i].wedge)
          {
            a.push_back(make_pair(i,int(avg_y)));
            n++;
          }
      }
  std::sort(a.begin(), a.end(),comp_thickness);
  double t;
  if (n > 0)
    t = a[(n - 1) / 2].second;
  else
    t = 1.5;
  //cout << "----------------" << endl;
  if (!triangle)
    for (int i=0; i<a.size(); i++)
      if (a[i].second > 2*t && a[i].first >= 0)
	bond[a[i].first].wedge = true;
  return (t);
}

void collapse_double_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double dist)
{
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && bond[i].type == 2 && bond[i].conjoined) // uninitialized value here!!!
      for (int j = 0; j < n_bond; j++)
        if (bond[j].exists && j != i && bond[j].type == 1 && bond_length(bond, j, atom) <= dist)
          {
            if (bond[j].a == bond[i].a)
              {
                bond[j].exists = false;
                atom[bond[i].a].x = (atom[bond[i].a].x + atom[bond[j].b].x) / 2;
                atom[bond[i].a].y = (atom[bond[i].a].y + atom[bond[j].b].y) / 2;

		atom[bond[i].a].min_x = min(atom[bond[i].a].min_x,atom[bond[j].b].min_x);
                atom[bond[i].a].min_y = min(atom[bond[i].a].min_y,atom[bond[j].b].min_y);
		atom[bond[i].a].max_x = max(atom[bond[i].a].max_x,atom[bond[j].b].max_x);
                atom[bond[i].a].max_y = max(atom[bond[i].a].max_y,atom[bond[j].b].max_y);
                for (int k = 0; k < n_bond; k++)
                  if (bond[k].exists)
                    {
                      if (bond[k].a == bond[j].b)
                        {
                          bond[k].a = bond[i].a;
                        }
                      else if (bond[k].b == bond[j].b)
                        {
                          bond[k].b = bond[i].a;
                        }
                    }
              }
            else if (bond[j].b == bond[i].a)
              {
                bond[j].exists = false;
                atom[bond[i].a].x = (atom[bond[i].a].x + atom[bond[j].a].x) / 2;
                atom[bond[i].a].y = (atom[bond[i].a].y + atom[bond[j].a].y) / 2;

		atom[bond[i].a].min_x = min(atom[bond[i].a].min_x,atom[bond[j].a].min_x);
                atom[bond[i].a].min_y = min(atom[bond[i].a].min_y,atom[bond[j].a].min_y);
		atom[bond[i].a].max_x = max(atom[bond[i].a].max_x,atom[bond[j].a].max_x);
                atom[bond[i].a].max_y = max(atom[bond[i].a].max_y,atom[bond[j].a].max_y);

                for (int k = 0; k < n_bond; k++)
                  if (bond[k].exists)
                    {
                      if (bond[k].a == bond[j].a)
                        {
                          bond[k].a = bond[i].a;
                        }
                      else if (bond[k].b == bond[j].a)
                        {
                          bond[k].b = bond[i].a;
                        }
                    }
              }
            else if (bond[j].a == bond[i].b)
              {
                bond[j].exists = false;
                atom[bond[i].b].x = (atom[bond[i].b].x + atom[bond[j].b].x) / 2;
                atom[bond[i].b].y = (atom[bond[i].b].y + atom[bond[j].b].y) / 2;
		atom[bond[i].b].min_x = min(atom[bond[i].b].min_x,atom[bond[j].b].min_x);
                atom[bond[i].b].min_y = min(atom[bond[i].b].min_y,atom[bond[j].b].min_y);
		atom[bond[i].b].max_x = max(atom[bond[i].b].max_x,atom[bond[j].b].max_x);
                atom[bond[i].b].max_y = max(atom[bond[i].b].max_y,atom[bond[j].b].max_y);
                for (int k = 0; k < n_bond; k++)
                  if (bond[k].exists)
                    {
                      if (bond[k].a == bond[j].b)
                        {
                          bond[k].a = bond[i].b;
                        }
                      else if (bond[k].b == bond[j].b)
                        {
                          bond[k].b = bond[i].b;
                        }
                    }
              }
            else if (bond[j].b == bond[i].b)
              {
                bond[j].exists = false;
                atom[bond[i].b].x = (atom[bond[i].b].x + atom[bond[j].a].x) / 2;
                atom[bond[i].b].y = (atom[bond[i].b].y + atom[bond[j].a].y) / 2;
		atom[bond[i].b].min_x = min(atom[bond[i].b].min_x,atom[bond[j].a].min_x);
                atom[bond[i].b].min_y = min(atom[bond[i].b].min_y,atom[bond[j].a].min_y);
		atom[bond[i].b].max_x = max(atom[bond[i].b].max_x,atom[bond[j].a].max_x);
                atom[bond[i].b].max_y = max(atom[bond[i].b].max_y,atom[bond[j].a].max_y);
                for (int k = 0; k < n_bond; k++)
                  if (bond[k].exists)
                    {
                      if (bond[k].a == bond[j].a)
                        {
                          bond[k].a = bond[i].b;
                        }
                      else if (bond[k].b == bond[j].a)
                        {
                          bond[k].b = bond[i].b;
                        }
                    }
              }
          }
}

void find_up_down_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double thickness)
{
  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists && bond[i].type == 2)
      {
        if (atom[bond[i].a].x > atom[bond[i].b].x)
          bond_end_swap(bond, i);
        if (atom[bond[i].a].x == atom[bond[i].b].x && atom[bond[i].a].y > atom[bond[i].b].y)
          bond_end_swap(bond, i);

        for (int j = 0; j < n_bond; j++)
          if (bond[j].exists && bond[j].type == 1 && !bond[j].wedge && !bond[j].hash)
            {
              bond[j].down = false;
              bond[j].up = false;
              if (bond[j].b == bond[i].a)
                {
                  double h = distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                  atom[bond[i].b].y, atom[bond[j].a].x, atom[bond[j].a].y);
                  if (h > thickness)
                    bond[j].down = true;
                  else if (h < -thickness)
                    bond[j].up = true;
                }
              else if (bond[j].a == bond[i].a)
                {
                  bond_end_swap(bond, j);
                  double h = distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                  atom[bond[i].b].y, atom[bond[j].a].x, atom[bond[j].a].y);
                  if (h > thickness)
                    bond[j].down = true;
                  else if (h < -thickness)
                    bond[j].up = true;
                }
              else if (bond[j].a == bond[i].b)
                {
                  double h = distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                  atom[bond[i].b].y, atom[bond[j].b].x, atom[bond[j].b].y);
                  if (h > thickness)
                    bond[j].up = true;
                  else if (h < -thickness)
                    bond[j].down = true;
                }
              else if (bond[j].b == bond[i].b)
                {
                  bond_end_swap(bond, j);
                  double h = distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                  atom[bond[i].b].y, atom[bond[j].b].x, atom[bond[j].b].y);
                  if (h > thickness)
                    bond[j].up = true;
                  else if (h < -thickness)
                    bond[j].down = true;
                }
            }
      }
}


void find_old_aromatic_bonds(const potrace_path_t *p, vector<bond_t> &bond, int n_bond, vector<atom_t> &atom,
                             int n_atom, double avg)
{
  const potrace_path_t *p1 = p;

  for (int i = 0; i < n_bond; i++)
    if (bond[i].exists)
      bond[i].arom = false;
  while (p != NULL)
    {
      if ((p->sign == int('-')) && detect_curve(bond, n_bond, p))
        {
          potrace_path_t *child = p->childlist;
          if (child != NULL && child->sign == int('+'))
            {
              potrace_path_t *gchild = child->childlist;
              if (gchild != NULL && gchild->sign == int('-'))
                {
                  for (int i = 0; i < n_bond; i++)
                    if (bond[i].exists && bond[i].curve == p)
                      bond[i].arom = true;
                  delete_curve_with_children(atom, bond, n_atom, n_bond, child);
                }
            }
        }
      p = p->next;
    }

  while (p1 != NULL)
    {
      if (p1->sign == int('+') && detect_curve(bond, n_bond, p1))
        {
          potrace_path_t *child = p1->childlist;
          if (child != NULL && child->sign == int('-'))
            {
              vector<int> vert;
              double circum = 0;
              for (int i = 0; i < n_bond; i++)
                if (bond[i].exists && bond[i].curve == p1)
                  circum += bond_length(bond, i, atom);
              for (int i = 0; i < n_atom; i++)
                if (atom[i].exists && atom[i].curve == p1)
                  vert.push_back(i);
              if (vert.size() > 4)
                {
                  double diameter = 0, center_x = 0, center_y = 0;
                  int num = 0;
                  for (unsigned int i = 0; i < vert.size(); i++)
                    {
                      for (unsigned int j = i + 1; j < vert.size(); j++)
                        {
                          double dist = distance(atom[vert[i]].x, atom[vert[i]].y, atom[vert[j]].x, atom[vert[j]].y);
                          if (dist > diameter)
                            diameter = dist;
                        }
                      center_x += atom[vert[i]].x;
                      center_y += atom[vert[i]].y;
                      num++;
                    }
                  center_x /= num;
                  center_y /= num;
                  bool centered = true;
                  for (unsigned int i = 0; i < vert.size(); i++)
                    {
                      double dist = distance(atom[vert[i]].x, atom[vert[i]].y, center_x, center_y);
                      if (fabs(dist - diameter / 2) > V_DISPLACEMENT)
                        centered = false;
                    }

                  if (circum < PI * diameter && diameter > avg / 2 && diameter < 3 * avg && centered)
                    {
                      delete_curve_with_children(atom, bond, n_atom, n_bond, p1);
                      for (int i = 0; i < n_bond; i++)
                        if (bond[i].exists)
                          {
                            double dist = distance((atom[bond[i].a].x + atom[bond[i].b].x) / 2, (atom[bond[i].a].y
                                                   + atom[bond[i].b].y) / 2, center_x, center_y);
                            double ang = angle4(atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[i].a].x,
                                                atom[bond[i].a].y, center_x, center_y, atom[bond[i].a].x, atom[bond[i].a].y);
                            ang = acos(ang) * 180.0 / PI;
                            if (ang < 90 && dist < (avg / 3 + diameter / 2))
                              {
                                bond[i].arom = true;
                              }
                          }
                    }

                }
            }
        }
      p1 = p1->next;
    }

}

void flatten_bonds(vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, double maxh)
{
  bool found = true;

  while (found)
    {
      found = false;
      for (int i = 0; i < n_bond; i++)
        if (bond[i].exists && bond[i].type < 3)
          {
            int n = 0;
            int f = i;
            double li = bond_length(bond, i, atom);
            if (atom[bond[i].a].label == " ")
              {
                for (int j = 0; j < n_bond; j++)
                  if (j != i && bond[j].exists && bond[j].type < 3 && (bond[i].a == bond[j].a || bond[i].a
                      == bond[j].b))
                    {
                      n++;
                      f = j;
                    }
                double lf = bond_length(bond, f, atom);
                if (n == 1)
                  {
                    if (bond[i].a == bond[f].b)
                      {
                        double h = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y,
                                                             atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[f].a].x, atom[bond[f].a].y));
                        double d = distance_from_bond_x_a(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                          atom[bond[i].b].y, atom[bond[f].a].x, atom[bond[f].a].y);
                        if (h <= maxh && d < 0)
                          {
                            bond[f].exists = false;
                            atom[bond[f].b].exists = false;
                            bond[i].a = bond[f].a;
                            if (lf > li)
                              bond[i].type = bond[f].type;
                            if (bond[f].arom)
                              bond[i].arom = true;
                            if (bond[f].hash)
                              bond[i].hash = true;
                            if (bond[f].wedge)
                              bond[i].wedge = true;
                            found = true;
                          }
                      }
                    else
                      {
                        double h = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y,
                                                             atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[f].b].x, atom[bond[f].b].y));
                        double d = distance_from_bond_x_a(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                          atom[bond[i].b].y, atom[bond[f].b].x, atom[bond[f].b].y);
                        if (h <= maxh && d < 0)
                          {
                            bond[f].exists = false;
                            atom[bond[f].a].exists = false;

                            if (bond[f].hash || bond[f].wedge)
                              {
                                bond[i].a = bond[i].b;
                                bond[i].b = bond[f].b;
                              }
                            else
                              bond[i].a = bond[f].b;
                            if (lf > li)
                              bond[i].type = bond[f].type;
                            if (bond[f].arom)
                              bond[i].arom = true;
                            if (bond[f].hash)
                              bond[i].hash = true;
                            if (bond[f].wedge)
                              bond[i].wedge = true;
                            found = true;
                          }
                      }
                  }
              }

            n = 0;
            f = i;
            if (atom[bond[i].b].label == " ")
              {
                for (int j = 0; j < n_bond; j++)
                  if (j != i && bond[j].exists && bond[j].type < 3 && (bond[i].b == bond[j].a || bond[i].b == bond[j].b))
                    {
                      n++;
                      f = j;
                    }
                double lf = bond_length(bond, f, atom);
                if (n == 1)
                  {
                    if (bond[i].b == bond[f].b)
                      {
                        double h = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y,
                                                             atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[f].a].x, atom[bond[f].a].y));
                        double d = distance_from_bond_x_b(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                          atom[bond[i].b].y, atom[bond[f].a].x, atom[bond[f].a].y);
                        if (h <= maxh && d > 0)
                          {
                            bond[f].exists = false;
                            atom[bond[f].b].exists = false;
                            if (bond[f].hash || bond[f].wedge)
                              {
                                bond[i].b = bond[i].a;
                                bond[i].a = bond[f].a;
                              }
                            else
                              bond[i].b = bond[f].a;
                            if (lf > li)
                              bond[i].type = bond[f].type;
                            if (bond[f].arom)
                              bond[i].arom = true;
                            if (bond[f].hash)
                              bond[i].hash = true;
                            if (bond[f].wedge)
                              bond[i].wedge = true;
                            found = true;
                          }
                      }
                    else
                      {
                        double h = fabs(distance_from_bond_y(atom[bond[i].a].x, atom[bond[i].a].y,
                                                             atom[bond[i].b].x, atom[bond[i].b].y, atom[bond[f].b].x, atom[bond[f].b].y));
                        double d = distance_from_bond_x_b(atom[bond[i].a].x, atom[bond[i].a].y, atom[bond[i].b].x,
                                                          atom[bond[i].b].y, atom[bond[f].b].x, atom[bond[f].b].y);
                        if (h <= maxh && d > 0)
                          {
                            bond[f].exists = false;
                            atom[bond[f].a].exists = false;
                            bond[i].b = bond[f].b;
                            if (lf > li)
                              bond[i].type = bond[f].type;
                            if (bond[f].arom)
                              bond[i].arom = true;
                            if (bond[f].hash)
                              bond[i].hash = true;
                            if (bond[f].wedge)
                              bond[i].wedge = true;
                            found = true;
                          }
                      }
                  }
              }
          }
    }
}


void mark_terminal_atoms(const vector<bond_t> &bond, int n_bond, vector<atom_t> &atom, int n_atom)
{
  for (int i = 0; i < n_atom; i++)
    atom[i].terminal = false;

  for (int j = 0; j < n_bond; j++)
    if (bond[j].exists && bond[j].type == 1 && !bond[j].arom)
      {
        if (terminal_bond(bond[j].a, j, bond, n_bond))
          atom[bond[j].a].terminal = true;
        if (terminal_bond(bond[j].b, j, bond, n_bond))
          atom[bond[j].b].terminal = true;
      }
}



void find_limits_on_avg_bond(double &best_bond, const vector<vector<double> > &pages_of_avg_bonds,
                             const vector<vector<double> > &pages_of_ind_conf)
{
  double max_ind_conf = -FLT_MAX;
  
  for (unsigned int l = 0; l < pages_of_ind_conf.size(); l++)
    for (unsigned int i = 0; i < pages_of_ind_conf[l].size(); i++)
      if (max_ind_conf < pages_of_ind_conf[l][i])
        {
          max_ind_conf = pages_of_ind_conf[l][i];
          best_bond = pages_of_avg_bonds[l][i];
        }
}

void remove_small_bonds_in_chars(vector<atom_t> &atom, vector<bond_t> &bond,vector<letters_t> &letters)
{
  for (int i = 0; i < bond.size(); i++)
    if (bond[i].exists)
      for (int j=0; j < letters.size(); j++)
	if (atom[bond[i].a].x >= letters[j].min_x && atom[bond[i].a].x <= letters[j].max_x &&
	    atom[bond[i].a].y >= letters[j].min_y && atom[bond[i].a].y <= letters[j].max_y &&
	    atom[bond[i].b].x >= letters[j].min_x && atom[bond[i].b].x <= letters[j].max_x &&
	    atom[bond[i].b].y >= letters[j].min_y && atom[bond[i].b].y <= letters[j].max_y)
	  bond[i].exists = false;
}
