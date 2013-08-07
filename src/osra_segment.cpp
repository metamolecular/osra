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
// File: osra_segment.cpp
//
// Defines page segmentation functions
//

#include <iostream> // std::ostream, std::cout

#include "osra.h"
#include "osra_common.h"
#include "osra_segment.h"
#include "osra_labels.h"
#include "osra_ocr.h"

unsigned int distance_between_points(const point_t &p1, const point_t &p2)
{
  return max(abs(p1.x - p2.x), abs(p1.y - p2.y));
}

unsigned int distance_between_segments(const vector<point_t> &s1, const vector<point_t> &s2)
{
  int r = INT_MAX;
  int d;

  for (vector<point_t>::const_iterator i = s1.begin(); i != s1.end(); i++)
    for (vector<point_t>::const_iterator j = s2.begin(); j != s2.end(); j++)
      {
        d = distance_between_points(*i, *j);
        if (d < r)
          r = d;
      }

  /*
  unsigned int ii, jj;
  #pragma omp parallel
  {
  	int priv_min = INT_MAX;
  #pragma omp for
  	for (ii = 0; ii < s1.size(); ii++) {
  		for (jj = 0; jj < s2.size(); jj++) {
  			int d = distance_between_points(s1[ii], s2[jj]);
  			if (d < priv_min)
  				priv_min = d;
  		}
  	}

  	//#pragma omp flush (r)
  	if (priv_min < r) {
  #pragma omp critical
  		{
  			if (priv_min < r)
  				r = priv_min;
  		}
  	}
  }
  */

  return r;
}

void find_connected_components(const Image &image, double threshold, const ColorGray &bgColor,
                               vector<list<point_t> > &segments, vector<vector<point_t> > &margins, bool adaptive)
{
  point_t p;
  list<point_t> points;
  int speckle_area = 2;
  if (adaptive)
    {
      int speckle_side = min(image.columns(), image.rows()) / 200;
      speckle_area = speckle_side * speckle_side;
      if (speckle_area < 2) speckle_area = 2;
    }

  vector<vector<int> > tmp(image.columns(), vector<int> (image.rows(), 0));

  for (unsigned int i = 0; i < image.columns(); i++)
    for (unsigned int j = 0; j < image.rows(); j++)
      if (get_pixel(image, bgColor, i, j, threshold) == 1) // populate with low threshold for future anisotropic smoothing
        tmp[i][j] = 1;


  for (unsigned int i = 0; i < image.columns(); i++)
    for (unsigned int j = 0; j < image.rows(); j++)
      if (tmp[i][j] == 1)
        {
          tmp[i][j] = 2;
          p.x = i;
          p.y = j;
          points.push_back(p);
          list<point_t> new_segment;
          vector<point_t> new_margin;
          int counter = 0;
          point_t p1;
          while (!points.empty())
            {
              p = points.back();
              points.pop_back();
              new_segment.push_back(p);
              tmp[p.x][p.y] = -1;
              bool on_the_margin = false;

              // "k" should be in range "[0 .. image.columns) intercepted with [p.x - 1 .. p.x + 2)" ==> "p.x + 2" should be positive ==> "p.x >= -1"
              // "l" should be in range "[0 .. image.rows)    intercepted with [p.y - 1 .. p.y + 2)" ==> "p.y + 2" should be positive ==> "p.y >= -1"
              if (p.x >= -1 && p.y >= -1)
                {
                  unsigned int x_lower = p.x > 1 ? p.x - 1 : 0; // "k" cannot be less then zero
                  unsigned int y_lower = p.y > 1 ? p.y - 1 : 0; // "l" cannot be less then zero
                  unsigned int x_upper = p.x + 2;
                  if (x_upper > image.columns())
                    x_upper = image.columns();
                  unsigned int y_upper = p.y + 2;
                  if (y_upper > image.rows())
                    y_upper = image.rows();

                  for (int k = x_lower; k < x_upper; k++)
                    for (int l = y_lower; l < y_upper; l++)
                      {
                        if (tmp[k][l] == 1)
                          {
                            p1.x = k;
                            p1.y = l;
                            points.push_back(p1);
                            tmp[k][l] = 2;
                          }
                        else if ((int) k != p.x && (int) l != p.y && tmp[k][l] == 0)
                          {
                            on_the_margin = true;
                          }
                      }
                }

              if (on_the_margin && (new_margin.size() < PARTS_IN_MARGIN || (counter % PARTS_IN_MARGIN) == 0))
                new_margin.push_back(p);
              if (on_the_margin)
                counter++;
            }
          if (segments.size() > MAX_SEGMENTS)
            return;
          if (new_segment.size() > speckle_area)
            {
              segments.push_back(new_segment);
              margins.push_back(new_margin);
            }
        }
}

unsigned int area_ratio(unsigned int a, unsigned int b)
{
  double r = max(a, b) / min(a, b);
  return (unsigned int) r;
}

void build_distance_matrix(const vector<vector<point_t> > &margins, unsigned int max_dist,
                           vector<vector<int> > &distance_matrix, vector<vector<int> > &features, const vector<list<point_t> > &segments,
                           unsigned int max_area_ratio, vector<vector<int> > &area_matrix)
{
  unsigned int d;
  unsigned int ar;

  for (unsigned int s1 = 0; s1 < margins.size(); s1++)
    for (unsigned int s2 = s1 + 1; s2 < margins.size(); s2++)
      if (distance_between_points(margins[s1].front(), margins[s2].front()) < (PARTS_IN_MARGIN
          * margins[s1].size() + PARTS_IN_MARGIN * margins[s2].size()) / 2 + max_dist)
        {
          d = distance_between_segments(margins[s1], margins[s2]);
          if (d < max_dist)
            {
              distance_matrix[s1][s2] = d;
              distance_matrix[s2][s1] = d;
              ar = area_ratio(segments[s1].size(), segments[s2].size());
              //cout << ar << endl;
              area_matrix[s1][s2] = ar;
              area_matrix[s2][s1] = ar;
              if (ar < max_area_ratio && d < max_dist)
                features[ar][d]++;
            }
        }
}

list<list<list<point_t> > > build_explicit_clusters(const list<list<int> > &clusters,
    const vector<list<point_t> > &segments)
{
  list<list<list<point_t> > > explicit_clusters;
  for (list<list<int> >::const_iterator c = clusters.begin(); c != clusters.end(); c++)
    {
      list<list<point_t> > set_of_segments;
      for (list<int>::const_iterator s = c->begin(); s != c->end(); s++)
        if (!segments[*s].empty())
          set_of_segments.push_back(segments[*s]);
      if (!set_of_segments.empty())
        explicit_clusters.push_back(set_of_segments);
    }
  return explicit_clusters;
}

void remove_separators(vector<list<point_t> > &segments, vector<vector<point_t> > &margins, double max_aspect,
                       unsigned int size)
{
  vector<list<point_t> >::iterator s;
  vector<vector<point_t> >::iterator m;
  s = segments.begin();
  m = margins.begin();

  while (s != segments.end() && m != margins.end())
    {
      if (s->size() <= size)
        {
          s++;
          m++;
          continue;
        }

      int stop = INT_MAX, sleft = INT_MAX, sbottom = 0, sright = 0;
      for (list<point_t>::iterator p = s->begin(); p != s->end(); p++)
        {
          if (p->x < sleft)
            sleft = p->x;
          if (p->x > sright)
            sright = p->x;
          if (p->y < stop)
            stop = p->y;
          if (p->y > sbottom)
            sbottom = p->y;
        }
      double aspect = 0;

      if (sright != sleft)
        aspect = 1. * (sbottom - stop+1) / (sright - sleft+1); // where did right and left come from?
      if (aspect > max_aspect || aspect < 1. / max_aspect)
        {
          s = segments.erase(s);
          m = margins.erase(m);
        }
      else
        {
          s++;
          m++;
        }
    }
}

void remove_tables_old(vector<list<point_t> > &segments, vector<vector<point_t> > &margins, unsigned int size)
{
  vector<list<point_t> >::iterator s;
  vector<vector<point_t> >::iterator m;
  s = segments.begin();
  m = margins.begin();

  while (s != segments.end() && m != margins.end())
    {
      if (m->size() <= size)
        {
          s++;
          m++;
          continue;
        }

      int top = INT_MAX, left = INT_MAX, bottom = 0, right = 0;
      int border_count = 0;
      for (vector<point_t>::iterator p = m->begin(); p != m->end(); p++)
        {
          if (p->x < left)
            left = p->x;
          if (p->x > right)
            right = p->x;
          if (p->y < top)
            top = p->y;
          if (p->y > bottom)
            bottom = p->y;
        }

      double aspect = FLT_MAX;
      if (right != left)
        aspect = 1. * (bottom - top) / (right - left);
      if (aspect >= MAX_ASPECT || aspect <= MIN_ASPECT)
        {
          s++;
          m++;
          continue;
        }

      for (vector<point_t>::iterator p = m->begin(); p != m->end(); p++)
        if (p->x - left < 2 || right - p->x < 2 || p->y - top < 2 || bottom - p->y < 2)
          {
            border_count++;
          }

      //cout << border_count << " " << 2 * (right - left) + 2 * (bottom - top) << endl;

      if (border_count > BORDER_COUNT)
        {
          s = segments.erase(s);
          m = margins.erase(m);
        }
      else
        {
          s++;
          m++;
        }
    }
}

double cos_angle_between_points(point_t a, point_t b, point_t c)
{
  double v1_x = a.x-b.x;
  double v1_y = a.y-b.y;
  double v2_x = c.x-b.x;
  double v2_y = c.y-b.y;
  double v1 = sqrt(v1_x*v1_x+v1_y*v1_y);
  double v2 = sqrt(v2_x*v2_x+v2_y*v2_y);
  if (v1>0 && v2>0)
    return (v1_x*v2_x+v1_y*v2_y)/(v1*v2);
  else
    return FLT_MAX;
}

// clockwise actual angle is considered positive
pair<double,double> find_rotation(point_t top, point_t left, point_t bottom, point_t right)
{
  double top_angle = fabs(cos_angle_between_points(left,top,right));
  double right_angle = fabs(cos_angle_between_points(top,right,bottom));
  double bottom_angle = fabs(cos_angle_between_points(right,bottom,left));
  double left_angle = fabs(cos_angle_between_points(bottom,left,top));
  double dx,dy;
  if (top_angle<right_angle && top_angle<left_angle && top_angle<bottom_angle)
    {
      // angle at the top point is the closest to 90 degrees
      if (top.x-left.x < right.x - top.x)
	{
	  // rotation is clockwise
	  dx = top.x-left.x;
	  dy = left.y-top.y;
	}
      else
	{
	  // rotation is counter-clockwise
	  dx = top.x-right.x;
	  dy = right.y-top.y;
	}
    }
  else if (left_angle<right_angle && left_angle<top_angle && left_angle<bottom_angle)
    {
      // angle at the left is the closest to 90 degrees 
      if (bottom.y-left.y < left.y - top.y)
	{
	  // clockwise
	  dx = top.x-left.x;
	  dy = left.y-top.y;
	} 
      else
	{
	  // counter-clockwise
	  dx = left.x-bottom.x;
	  dy = bottom.y-left.y;
	}
    }
  else if (right_angle<top_angle && right_angle<left_angle && right_angle<bottom_angle)
    {
      // angle at the right is the closest to 90 degrees
      if (right.y - top.y < bottom.y - right.y)
	{
	  // clockwise
	  dx = right.x-bottom.x;
	  dy = bottom.y-right.y;
	}
      else
	{
	  // counter-clockwise
	  dx = top.x-right.x;
	  dy = right.y-top.y;
	}
    }
  else
    {
      // assume that angle at the bottom is the closest to 90 degrees
      if (right.x-bottom.x < bottom.x-left.x)
	{
	  // clockwise
	  dx = right.x-bottom.x;
	  dy = bottom.y-right.y;
	}
      else 
	{
	  // counter-clockwise
	  dx = left.x - bottom.x;
	  dy = bottom.y - left.y;
	}
    }
  double s = dx / sqrt(dx*dx+dy*dy);
  double c = dy / sqrt(dx*dx+dy*dy);
  return(make_pair(s,c));
}

int border_count_in_rotated_frame(vector<vector<point_t> >::iterator m,
				  point_t top_point,point_t left_point,point_t bottom_point,point_t right_point,
				  pair <double,double> &sin_cos)
{
  double s = -sin_cos.first;
  double c = sin_cos.second;

  double x1 = top_point.x*c-top_point.y*s;
  double y1 = top_point.x*s+top_point.y*c;
  double x2 = left_point.x*c-left_point.y*s;
  double y2 = left_point.x*s+left_point.y*c;
  double x3 = bottom_point.x*c-bottom_point.y*s;
  double y3 = bottom_point.x*s+bottom_point.y*c;
  double x4 = right_point.x*c-right_point.y*s;
  double y4 = right_point.x*s+right_point.y*c;

  double left = min(min(x1,x2),min(x3,x4));
  double right = max(max(x1,x2),max(x3,x4));
  double top = min(min(y1,y2),min(y3,y4));
  double bottom = max(max(y1,y2),max(y3,y4));

  int border_count = 0;
  for (vector<point_t>::iterator p = m->begin(); p != m->end(); p++)
    {
      double x = c*(p->x) - s*(p->y);
      double y = s*(p->x) + c*(p->y);
      
      if ((x - left)<2 || (right - x) < 2 || (y - top) < 2 || (bottom - y) < 2)
	border_count++;
    }
  return border_count;
}


void remove_tables(vector<list<point_t> > &segments, vector<vector<point_t> > &margins, unsigned int size)
{
  vector<list<point_t> >::iterator s;
  vector<vector<point_t> >::iterator m;
  s = segments.begin();
  m = margins.begin();

  while (s != segments.end() && m != margins.end())
    {
      if (m->size() <= size)
        {
          s++;
          m++;
          continue;
        }

      int top = INT_MAX, left = INT_MAX, bottom = 0, right = 0;
      point_t left_point,top_point,right_point,bottom_point;
      for (vector<point_t>::iterator p = m->begin(); p != m->end(); p++)
        {
          if (p->x < left)
	    {
	      left = p->x;
	      left_point = *p;
	    }
          if (p->x > right)
	    {
	      right = p->x;
	      right_point = *p;
	    }
          if (p->y < top)
	    {
	      top = p->y;
	      top_point = *p;
	    }
          if (p->y > bottom)
	    {
	      bottom = p->y;
	      bottom_point = *p;
	    }
        }

      double aspect = FLT_MAX;
      if (right != left)
        aspect = 1. * (bottom - top) / (right - left);
      if (aspect >= MAX_ASPECT || aspect <= 1./MAX_ASPECT)
        {
          s++;
          m++;
          continue;
        }

      double area = s->size();
      double square_area = (bottom - top+1) * (right - left+1);
      double ratio = 0;
      if (square_area != 0)
	ratio = area / square_area;
      if (ratio > MAX_RATIO || ratio == 0)
	{
          s++;
          m++;
          continue;
        }
      pair<double,double> zero_angle = make_pair(0,1);
      int border_count = border_count_in_rotated_frame(m,top_point,left_point,bottom_point,right_point,zero_angle);

      if (PARTS_IN_MARGIN*border_count > BORDER_COUNT)
        {
          s = segments.erase(s);
          m = margins.erase(m);
        }
      else
        {
	  // perform rotation
	  pair<double,double> sin_cos = find_rotation(top_point,left_point,bottom_point,right_point);
	  int rotated_border_count = border_count_in_rotated_frame(m,top_point,left_point,bottom_point,right_point,sin_cos);

	    if (PARTS_IN_MARGIN*rotated_border_count > BORDER_COUNT && fabs(sin_cos.first)<sin(10*PI/180) && sin_cos.second>cos(10*PI/180))
	    {
	      s = segments.erase(s);
	      m = margins.erase(m);
	    }
	  else
	    { 
	      s++;
	      m++;
	    }
        }
    }
}

list<list<int> > assemble_clusters(const vector<vector<point_t> > &margins, int dist,
                                   const vector<vector<int> > &distance_matrix, vector<int> &avail, bool text,
                                   const vector<vector<int> > &area_matrix)
{
  list<list<int> > clusters;
  list<int> bag;

  for (unsigned int s = 0; s < margins.size(); s++)
    if (avail[s] == 1)
      {
        bag.push_back(s);
        avail[s] = 2;
        list<int> new_cluster;
        while (!bag.empty())
          {
            int c = bag.back();
            bag.pop_back();
            new_cluster.push_back(c);
            avail[c] = 0;
            for (unsigned int i = 0; i < margins.size(); i++)
              if (avail[i] == 1 && distance_matrix[c][i] < dist)
                // && (!text || area_matrix[i][c] <= 10))
                {
                  bag.push_back(i);
                  avail[i] = 2;
                }
          }
        clusters.push_back(new_cluster);
      }

  return (clusters);
}

void remove_text_blocks(const list<list<int> > &clusters, const vector<list<point_t> > &segments, vector<int> &avail)
{
  for (list<list<int> >::const_iterator c = clusters.begin(); c != clusters.end(); c++)
    {
      unsigned int area = 0, square_area = 0;
      double ratio = 0, aspect = 0;
      int top = INT_MAX, left = INT_MAX, bottom = 0, right = 0;
      bool fill_below_max = false;

      for (list<int>::const_iterator i = c->begin(); i != c->end(); i++)
        if (!segments[*i].empty())
          {
            int stop = INT_MAX, sleft = INT_MAX, sbottom = 0, sright = 0;
            for (list<point_t>::const_iterator p = segments[*i].begin(); p != segments[*i].end(); p++)
              {
                if (p->x < sleft)
                  sleft = p->x;
                if (p->x > sright)
                  sright = p->x;
                if (p->y < stop)
                  stop = p->y;
                if (p->y > sbottom)
                  sbottom = p->y;
              }

            area = segments[*i].size();
            square_area = (sbottom - stop+1) * (sright - sleft+1);

            if (square_area != 0)
              ratio = 1. * area / square_area;

            if (ratio < MAX_RATIO && ratio > 0)
              fill_below_max = true;

            if (sleft < left)
              left = sleft;
            if (sright > right)
              right = sright;
            if (stop < top)
              top = stop;
            if (sbottom > bottom)
              bottom = sbottom;
          }

      if (c->size() > TEXT_LINE_SIZE)
        {
          if (right != left)
            aspect = 1. * (bottom - top) / (right - left);
          if (aspect < MIN_ASPECT || aspect > MAX_ASPECT || !fill_below_max)
            for (list<int>::const_iterator i = c->begin(); i != c->end(); i++)
              avail[*i] = -1;
        }
    }
}

int locate_first_min(const vector<int> &stats)
{
  int peak = 1;

  for (unsigned int j = 3; j < stats.size(); j++)
    if (stats[j] > stats[j - 1] && stats[j] > stats[j + 1])
      {
        peak = j;
        break;
      }
  int dist = peak;
  for (unsigned int j = peak; j < stats.size(); j++)
    if (stats[j] < stats[j - 1] && stats[j] < stats[j + 1])
      {
        dist = j;
        break;
      }
  return (dist);
}

int locate_max_entropy(const vector<vector<int> > &features, unsigned int max_area_ratio, unsigned int max_dist,
                       vector<int> &stats)
{
  vector<double> entropy(max_area_ratio, 0);

  for (unsigned int i = 1; i < max_area_ratio; i++)
    {
      int count = 0;
      for (unsigned int j = 2; j < max_dist; j++)
        if (features[i][j] == 0)
          count++;
        else
          stats[j]++;

      if (count > 0)
        {
          double probability = 1. * count / (max_dist - 2);
          entropy[i] -= probability * log(probability);
        }
    }
  int start_b = 1;
  for (unsigned int i = 2; i < max_area_ratio; i++)
    {
      if (entropy[i] > entropy[start_b])
        start_b = i;
    }
  return (start_b);
}


bool bulge(const point_t tail, const point_t head, const list<point_t> & seg)
{
  bool r = false;
  vector<int> y(max(abs(head.x-tail.x),abs(head.y-tail.y))+1,0);
  int n=y.size();

  if (n<10) return false;

  bool horizontal = false;
  if (abs(head.x-tail.x)>abs(head.y-tail.y)) horizontal = true;

  for (list<point_t>::const_iterator p=seg.begin(); p!=seg.end(); p++)
    {
      int d;
      if (horizontal)
	d=abs(p->x-tail.x);
      else
	d=abs(p->y-tail.y);
      if (d<n) y[d]++;
    }

  int top=0;
  int pos=0;
  for (int i=0; i<n; i++)
    if (y[i]>top)
      {
	top = y[i];
	pos = i;
      }

  /*for (int i=0; i<n; i++)
   cout<<y[i]<<" ";
  cout<<endl;
  */
  if (pos<3) return false;
  int midpoint = min(int(0.75*n),pos-3);
  if (midpoint<0) return false;
 
  double avg=0;
  for (int i=0; i<midpoint; i++)
    avg +=y[i];
  avg /=int(midpoint);

  bool flat = true;
  for (int i=0; i<midpoint; i++)
    if (fabs(y[i]-avg)>2) flat = false;

  bool left = true;
  for (int i=pos-1; i>=midpoint; i--)
    if (y[i]>y[i+1]+2) 
      left=false;
  bool right = true;
  for (int i=pos+1; i<n; i++)
    if (y[i]>y[i-1]+2) right=false;
  bool peak = true;
  if (top<1.5*avg || top-avg<2 || n-pos<3 || top-y[n-1]<2 || pos<n/2 || pos<5) peak = false;

  //cout<<flat<<" "<<left<<" "<<right<<" "<<peak<<" "<<pos<<endl;

  return flat && left && right && peak;
}


void find_arrows_pluses(vector<vector<point_t> > &margins, vector<list<point_t> > &segments, vector<arrow_t> &arrows, vector<plus_t> &pluses)
{
  const int len=50;
  for (int i=0; i<margins.size(); i++)
    {
      vector<int> hist(len,0);
      int top_pos=0;
      int top_value=0;
      point_t head, tail,center;
      int min_x, min_y, max_x, max_y;

      if (segments[i].size()>1000)
	build_hist(margins[i],hist,len,top_pos,top_value,head,tail,center,min_x, min_y, max_x, max_y);
      else
	build_hist(segments[i],hist,len,top_pos,top_value,head,tail,center,min_x, min_y, max_x, max_y);

      if (top_value>5)
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
	      bool ba=bulge(tail,head,segments[i]);
	      bool bb=bulge(head,tail,segments[i]);
	      //cout<<tail.x<<" "<<tail.y<<" "<<ba<<" "<<bb<<endl;
	      if (ba || bb)
		{
		  // we found an arrow!
		  arrow_t arrow(head,tail,min_x,min_y,max_x,max_y);
		  if (bb)
		    {
		      arrow.head = tail;
		      arrow.tail = head;
		    }
		  arrows.push_back(arrow);
		  margins[i].clear();
		  segments[i].clear();
		}
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
		  plus_t plus;
		  plus.center = center;
		  plus.min_x = min_x;
		  plus.min_y = min_y;
		  plus.max_x = max_x;
		  plus.max_y = max_y;
		  pluses.push_back(plus);
		}
	    }
	  
	}
    }
  vector<vector<point_t> >::iterator k = margins.begin();
  while (k!=margins.end())
    {
      if (k->empty()) k = margins.erase(k);
      else k++;
    }
 vector<list<point_t> >::iterator l = segments.begin();
  while (l!=segments.end())
    {
      if (l->empty()) l = segments.erase(l);
      else l++;
    }
  
}

bool comp_labels(const label_t &left, const label_t &right)
{
  if (left.x2 < right.x1)
    return (true);
  if (max(left.y1,left.y2) < min(right.y1,right.y2))
    return (true);
  return (false);
}

int comp_labels_int(const void *l, const void *r)
{
  label_t *left = (label_t *) l;
  label_t *right = (label_t *) r;
  if (left->x2 < right->x1)
    return (-1);
  if (max(left->y1,left->y2) < min(right->y1,right->y2))
    return (-1);
  return (1);
}

string ocr_agent_strings(const vector<list<point_t> >  &agents,const Image &image, double threshold, const ColorGray &bgColor, bool verbose)
{ 
  string agent_string;
  vector<letters_t> letters;
  for (int i=0; i<agents.size(); i++)
    {
      int left=INT_MAX;
      int right=0;
      int top=INT_MAX;
      int bottom=0;
      for (list<point_t>::const_iterator a=agents[i].begin(); a!=agents[i].end(); a++)
	{
	  if (a->x<left) left=a->x;
	  if (a->x>right) right=a->x;
	  if (a->y<top) top=a->y;
	  if (a->y>bottom) bottom=a->y;
	}
      if ((bottom - top) <= 2*MAX_FONT_HEIGHT && (right - left) <= 2*MAX_FONT_WIDTH && (bottom - top) > MIN_FONT_HEIGHT)
            {
              char label = 0;
              label = get_atom_label(image, bgColor, left, top, right, bottom, threshold, (right + left) / 2, top, true, verbose);
              if (label != 0)
                {
                  letters_t lt;
		  lt.a=label;
		  lt.x = (left + right) / 2;
		  lt.y  = (top + bottom) / 2;
		  lt.r  = distance(left, top, right, bottom) / 2;
		  lt.min_x = left;
		  lt.max_x = right;
		  lt.min_y = top;
		  lt.max_y = bottom;
		  lt.free = true;
                  letters.push_back(lt);
		}
	    }
    }
  vector<label_t> label;
  assemble_labels(letters, letters.size(), label);
 
  //sort(label.begin(),label.end(),comp_labels);

  if (!label.empty())
    qsort(&label[0],label.size(),sizeof(label_t),comp_labels_int);

  for (int i=0; i<label.size(); i++)
    if (!label[i].a.empty())
      agent_string += " "+label[i].a;
  return (agent_string);
}

void find_agent_strings(vector<vector<point_t> > &margins,vector<list<point_t> > &segments, vector<arrow_t> &arrows, 
			const Image &image, double threshold, const ColorGray &bgColor, bool verbose)
{
  for (int i=0; i<arrows.size(); i++)
    {
      vector<vector<point_t> >  agent_margins;
      vector<list<point_t> >  agents;

      double l=distance(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y);
      bool found=false;
      for (int j=0; j<margins.size(); j++)
	{
	  bool close=false;
	  bool within=true;
	  for (int k=0; k<margins[j].size(); k++)
	    { 
	      //     if (fabs(distance_from_bond_y(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,margins[j][k].x,margins[j][k].y))<MAX_FONT_HEIGHT) close=true;
	      //double d=distance_from_bond_x_a(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,margins[j][k].x,margins[j][k].y);
	      //if (d<0 || d>l) within=false;
	      if (fabs(distance((arrows[i].tail.x+arrows[i].head.x)/2,(arrows[i].tail.y+arrows[i].head.y)/2,margins[j][k].x,margins[j][k].y))<MAX_FONT_HEIGHT) close=true;
	    }
	  if (close && within)
	    {
	      agents.push_back(segments[j]);
	      agent_margins.push_back(margins[j]);
	      margins[j].clear();
	      segments[j].clear();
	      found=true;
	    }
	}

      while (found)
	{
	  found=false;

	  vector<vector<point_t> >::iterator k = margins.begin();
	  while (k!=margins.end())
	    {
	      if (k->empty()) k = margins.erase(k);
	      else k++;
	    }
	  vector<list<point_t> >::iterator l = segments.begin();
	  while (l!=segments.end())
	    {
	      if (l->empty()) l = segments.erase(l);
	      else l++;
	    }

	  for (int j=0; j<margins.size(); j++)
	    {
	      bool close=false;
	      for (int m=0; m<agent_margins.size(); m++)
		for (int k=0; k<margins[j].size(); k++)
		  for (int p=0; p<agent_margins[m].size(); p++)
		    if (distance(agent_margins[m][p].x,agent_margins[m][p].y,margins[j][k].x,margins[j][k].y)<MAX_FONT_HEIGHT/2) close=true;
	      if (close)
		{
		  agents.push_back(segments[j]);
		  agent_margins.push_back(margins[j]);
		  margins[j].clear();
		  segments[j].clear();
		  found=true;
		}
	    }
	}
      arrows[i].agent=ocr_agent_strings(agents,image,threshold,bgColor,verbose);
    }
}

list<list<list<point_t> > > find_segments(const Image &image, double threshold, const ColorGray &bgColor, bool adaptive, bool is_reaction, vector<arrow_t> &arrows, vector<plus_t> &pluses,
					  bool verbose)
{
  vector<list<point_t> > segments;
  vector<vector<point_t> > margins;
  list<list<list<point_t> > > explicit_clusters;

  // 1m34s

  find_connected_components(image, threshold, bgColor, segments, margins, adaptive);

  if (verbose)
    cout << "Number of segments: " << segments.size() << '.' << endl;

  if (segments.size() > MAX_SEGMENTS)
    {
      segments.clear();
      margins.clear();
    }
  if (is_reaction)
    {
      find_arrows_pluses(margins,segments,arrows, pluses);
      find_agent_strings(margins,segments,arrows,image,threshold,bgColor,verbose);
    }

  remove_separators(segments, margins, SEPARATOR_ASPECT, SEPARATOR_AREA);

  remove_tables(segments, margins, SEPARATOR_AREA);

  // 2m22s

  unsigned int max_dist = MAX_DIST;
  unsigned int max_area_ratio = MAX_AREA_RATIO;
  vector<vector<int> > distance_matrix(segments.size(), vector<int> (segments.size(), INT_MAX));
  vector<vector<int> > area_matrix(segments.size(), vector<int> (segments.size(), INT_MAX));
  vector<vector<int> > features(max_area_ratio, vector<int> (max_dist, 0));

  build_distance_matrix(margins, max_dist, distance_matrix, features, segments, max_area_ratio, area_matrix);

  // 2m53s

  vector<int> avail(margins.size(), 1);

  /*
  unsigned int ar;

  for (unsigned int i = 0; i < margins.size(); i++)
  	for (unsigned int j = i + 1; j < margins.size(); j++) {
  		ar = area_ratio(segments[i].size(), segments[j].size());
  		if (ar < max_area_ratio && distance_matrix[i][j] < max_dist)
  			features[ar][distance_matrix[i][j]]++;
  	}
  */

  // 5m53s -> new 4m15s

  vector<int> stats(max_dist, 0);
  int entropy_max = locate_max_entropy(features, max_area_ratio, max_dist, stats);

  int dist = SINGLE_IMAGE_DIST;

  if (entropy_max > THRESHOLD_LEVEL && !adaptive && margins.size() > 100)
    {
      vector<int> text_stats(max_dist, 0);
      for (unsigned int j = 2; j < max_dist; j++)
        {
          text_stats[j] = features[1][j];
          //cout << j << " " << text_stats[j] << endl;
        }

      int dist_text = locate_first_min(text_stats);

      const list<list<int> > &text_blocks = assemble_clusters(margins, dist_text, distance_matrix, avail, true,
                                            area_matrix);
      remove_text_blocks(text_blocks, segments, avail);

      dist = 2 * dist_text;
    }
  //  if (dist<20) dist = 20;
  for (unsigned int i = 0; i < margins.size(); i++)
    if (avail[i] != -1)
      avail[i] = 1;

  const list<list<int> > &clusters = assemble_clusters(margins, dist, distance_matrix, avail, false, area_matrix);

  explicit_clusters = build_explicit_clusters(clusters, segments);
  return explicit_clusters;
}

/*
void remove_brackets(int left, int right, int top, int bottom, list<list<list<point_t> > >::iterator c) {
	vector < vector<bool> > tmp(right - left + 1, vector<bool> (bottom - top + 1, false));
	vector < vector<bool> > global_pic(right - left + 1, vector<bool> (bottom - top + 1, false));

	for (list<list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
		for (list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
			global_pic[p->x - left][p->y - top] = true;

	bool found = false;
	//Image t(Geometry(right - left + 1, bottom - top + 1), "white");

	for (int i = left + FRAME; i < right - FRAME; i++)
		for (list<list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++) {
			vector<point_t> set;
			int x1 = INT_MAX, y1 = INT_MAX, x2 = 0, y2 = 0;
			for (list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
				if (p->x < i && i + (i - p->x) < right && global_pic[i + (i - p->x) - left][p->y - top]) {
					set.push_back(*p);
					if (p->x < x1)
						x1 = p->x;
					if (p->x > x2)
						x2 = p->x;
					if (p->y < y1)
						y1 = p->y;
					if (p->y > y2)
						y2 = p->y;
				}

			if (set.size() > 100 && (i - x2) > 5 && (x2 - x1) > 3 && (y2 - y1) > 5 && (x2 - x1) < (y2 - y1)) {
				int x = x2 - x1 + 1;
				int y = y2 - y1 + 1;
				int f = 1;
				if (y > 40)
					f = y / 40;
				x /= f;
				y /= f;

				unsigned char *pic = (unsigned char *) malloc(x * y);
				for (int j = 0; j < x * y; j++)
					pic[j] = 255;
				for (unsigned int p = 0; p < set.size(); p++)
					if ((set[p].y - y1) / f < y && (set[p].x - x1) / f < x)
						pic[((set[p].y - y1) / f) * x + (set[p].x - x1) / f] = 0;
				bool res = detect_bracket(x, y, pic);
				if (res) {
					//cout << set.size() << endl;
					for (unsigned int p = 0; p < set.size(); p++) {
						tmp[set[p].x - left][set[p].y - top] = true;
						tmp[i + (i - set[p].x) - left][set[p].y - top] = true;
						//t.pixelColor(set[p].x - left, set[p].y - top, "black");
						//t.pixelColor(i + (i - set[p].x) - left, set[p].y - top, "black");
					}
					found = true;
				}
			}
		}

	if (found) {
		//t.write("t.png");
		//exit(0);
		list<list<point_t> >::iterator s1 = c->begin();
		while (s1 != c->end()) {
			list<point_t>::iterator p1 = s1->begin();
			while (p1 != s1->end())
				if (tmp[p1->x - left][p1->y - top])
					p1 = s1->erase(p1);
				else
					p1++;
			if (s1->size() > 0)
				s1++;
			else
				s1 = c->erase(s1);
		}
	}
}
*/

int prune_clusters(list<list<list<point_t> > > &clusters, vector<box_t> &boxes)
{
  int n_boxes = 0;
  list<list<list<point_t> > >::iterator c = clusters.begin();

  while (c != clusters.end())
    {
      unsigned int area = 0, square_area = 0;
      double ratio = 0, aspect = 0;
      int top = INT_MAX, left = INT_MAX, bottom = 0, right = 0;
      bool fill_below_max = false;
      for (list<list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
        {
          int stop = INT_MAX, sleft = INT_MAX, sbottom = 0, sright = 0;
          for (list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
            {
              if (p->x < sleft)
                sleft = p->x;
              if (p->x > sright)
                sright = p->x;
              if (p->y < stop)
                stop = p->y;
              if (p->y > sbottom)
                sbottom = p->y;
            }

          area = s->size();
          square_area = (sbottom - stop+1) * (sright - sleft+1);
          ratio = 0;
          if (square_area != 0)
            ratio = 1. * area / square_area;

          if (ratio < MAX_RATIO && ratio > 0)
            fill_below_max = true;

          if (sleft < left)
            left = sleft;
          if (sright > right)
            right = sright;
          if (stop < top)
            top = stop;
          if (sbottom > bottom)
            bottom = sbottom;
        }

      if (right != left)
        aspect = 1. * (bottom - top) / (right - left);

      if (fill_below_max && aspect > MIN_ASPECT && aspect < MAX_ASPECT)
        {
          box_t b1;
          boxes.push_back(b1);
          boxes[n_boxes].x1 = left;
          boxes[n_boxes].y1 = top;
          boxes[n_boxes].x2 = right;
          boxes[n_boxes].y2 = bottom;

          //remove_brackets(left, right, top, bottom, c);

          for (list<list<point_t> >::const_iterator s = c->begin(); s != c->end(); s++)
            for (list<point_t>::const_iterator p = s->begin(); p != s->end(); p++)
              boxes[n_boxes].c.push_back(*p);
          c++;
          n_boxes++;
        }
      else
        {
          c = clusters.erase(c);
        }
    }
  return (n_boxes);
}

