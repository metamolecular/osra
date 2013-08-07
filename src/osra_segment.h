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
// Header: osra_segment.h
//
// Declares page segmentation functions
//

#ifndef OSRA_SEGMENT_H
#define OSRA_SEGMENT_H

#include <list> // sdt::list
#include <vector> // std::vector
#include <math.h> // fabs(double)
#include <float.h> // FLT_MAX
#include <limits.h> // INT_MAX
#include <algorithm> // std::min(double, double), std::max(double, double)
#include <Magick++.h>

using namespace std;
using namespace Magick;

// struct: point_s
//      a point of the image, used by image segmentation routines
struct point_s
{
  // int: x,y
  //    coordinates of the image point
  int x, y;
};
// typedef: point_t
//      defines point_t type based on point_s struct
typedef struct point_s point_t;

// struct: box_s
//      encompassing box structure for image segmentation
struct box_s
{
  // int: x1, y1, x2, y2
  //    coordinates of top-left and bottom-right corners
  int x1, y1, x2, y2;
  // array: c
  //    vector of points in the box
  vector<point_t> c;
};
// typedef: box_t
//      defines box_t type based on box_s struct
typedef struct box_s box_t;

// struct: arrow_s
// coordinates of tail and head of an arrow
struct arrow_s
{
arrow_s(point_t _head, point_t _tail,int _min_x,int _min_y,int _max_x,int _max_y) :
  head(_head),tail(_tail),min_x(_min_x),min_y(_min_y),max_x(_max_x),max_y(_max_y),linebreak(false),reversible(false),remove(false),agent("") {}
  arrow_s() {}
  // point_t: tail, head
  // tail and head of an arrow as points
  point_t tail,head;
  int min_x,min_y,max_x,max_y;
  string agent;
  bool linebreak;
  bool reversible;
  bool remove;
};
// typedef: arrow_t
// defines arrow_t type based on arrow_s struct
typedef struct arrow_s arrow_t;

struct plus_s
{
  point_t center;
   int min_x,min_y,max_x,max_y;
};
typedef struct plus_s plus_t;

//
// Section: Functions
//

// Function: find_segments()
//
// Performs page segmentation to different regions (text/graphics/linear etc.)
//
// Parameters:
// image - page image
// threshold - black-white binarization threshold
// bgColor - background color
// adaptive - flag set if adaptive thresholding has been used in grayscale conversion
// is_reaction - flag set if we're looking for reaction-specific symbols (arrows, plus signs etc.)
// arrows - a vector of arrows found during segmentation
// pluses - a vector of plus centers found during segmentation
// verbose - flag set for verbose reporting
//
// Returns:
// A list of clusters, each of which is a list of  connected segments each of which is a list of points
list<list<list<point_t> > > find_segments(const Image &image, double threshold, const ColorGray &bgColor, bool adaptive, bool is_reaction, vector<arrow_t> &arrows, vector<plus_t> &pluses, bool verbose);

// Function: prune_clusters()
//
// Prunes the list of clusters and retains only molecular structure images
//
// Parameters:
// clusters - a list of clusters detected by <find_segments()>
// boxes - a vector of <box_t> objects for molecular structure images
//
// Returns:
// Number of molecular structure images
int prune_clusters(list<list<list<point_t> > > &clusters, vector<box_t> &boxes);

template<class T>
void build_hist(const T &seg, vector<int> &hist, const int len, int &top_pos, int &top_value,point_t &head,point_t &tail, point_t &center, int &min_x, int &min_y, int &max_x, int &max_y)
{
  int l=seg.size();
  typename T::const_iterator j;
  center.x=0; center.y=0;
  min_x = INT_MAX;
  min_y = INT_MAX;
  max_x = 0;
  max_y = 0;
  for (j=seg.begin(); j!=seg.end(); j++)
    {
      center.x += j->x;
      center.y += j->y;
      min_x = min(min_x, j->x);
      min_y = min(min_y, j->y);
      max_x = max(max_x, j->x);
      max_y = max(max_y, j->y);
    }
  center.x /=l;  // Find the center of mass for the segment margin
  center.y /=l;

  for (j=seg.begin(); j!=seg.end(); j++)
    {
      int dx = j->x-center.x;
      int dy = j->y-center.y;
      double r=(double)sqrt(dx*dx+dy*dy);
      double theta=0.;
      if (dx!=0 || dy!=0)
	theta = atan2(dy,dx);
      int bin = (theta+M_PI)*len/(2*M_PI);
      if (bin>=len) bin -= len;
      hist[bin]++;          // build a histogram of occurencies in polar coordinates
      if (hist[bin]>=top_value)
	  {
	    top_pos = bin;                       // find the position of the highest peak
	    top_value = hist[bin];
	  }
    }

  double r_max=0;
  for (j=seg.begin(); j!=seg.end(); j++)
    {
      int dx = j->x-center.x;
      int dy = j->y-center.y;
      double r=(double)sqrt(dx*dx+dy*dy);
      double theta=0.;
      if (dx!=0 || dy!=0)
	theta = atan2(dy,dx);
      int bin = (theta+M_PI)*len/(2*M_PI);
      if (bin>=len) bin -= len;
      if (bin == top_pos && r>r_max)
	{
	  head = *j;
	  r_max = r;
	}
    }

  r_max=0;
  for (j=seg.begin(); j!=seg.end(); j++)
    {
      int dx = j->x-head.x;
      int dy = j->y-head.y;
      double r=(double)sqrt(dx*dx+dy*dy);
      if (r>r_max)
	{
	  r_max = r;
	  tail = *j;
	}
    }
}
#endif
