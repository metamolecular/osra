/******************************************************************************
 OSRA: Optical Structure Recognition

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

#include <sstream> // std:ostringstream
#include <string> // std:string
#include <vector> // std::vector
#include <set>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/reaction.h>

#include "osra_segment.h"
#include "osra_reaction.h"
#include "osra_common.h"

using namespace OpenBabel;
using namespace std;

//
// Create a reaction representation for input vector of structures
//
// Parameters:
//      page_of_structures - input vector of reactants, intermediates and products
//      output_format - format of the returned result, i.e. rsmi or cmlr
//
// Returns:
//      resulting reaction in the format set up by output_format parameter
//
string convert_page_to_reaction(const vector<string> &page_of_structures, const string &output_format, const vector <int> &reactants, const vector <int> &products, string value, bool reversible)
{
  string reaction;
  OBConversion *conv=new OBConversion;
  conv->SetInAndOutFormats(SUBSTITUTE_REACTION_FORMAT,output_format.c_str());
  ostringstream strstr;
  
  OBReaction react;
  for (int j=0; j<reactants.size(); j++)
    {
      shared_ptr<OBMol> reactant(new OBMol);
      conv->ReadString(reactant.get(), page_of_structures[reactants[j]]);
      react.AddReactant(reactant);
    }
  for (int j=0; j<products.size(); j++)
    {
      shared_ptr<OBMol> product(new OBMol);
      conv->ReadString(product.get(), page_of_structures[products[j]]);
      react.AddProduct(product);
    }
  //	  react.AddAgent(transition);
  if (reversible) react.SetReversible(true);

  trim(value);
  if (!value.empty())
    {
            OBPairData *label = new OBPairData;
            label->SetAttribute("OSRA_REACTION_AGENT");
            label->SetValue(value.c_str());
            react.SetData(label);
//      react.SetComment(value);
      react.SetTitle(value);
    }
  strstr << conv->WriteString(&react, true);
  reaction = strstr.str();
  if (output_format != "rxn") // rxn format seems to have a double-free problem in OB 2.3.1
   delete conv;
  return(reaction);
}

string convert_to_smiles_agent_structure(const string &structure)
{
  OBConversion conv;
  conv.SetInAndOutFormats(SUBSTITUTE_REACTION_FORMAT,"smi");
  string result;
  OBMol mol;
  if (conv.ReadString(&mol,structure))
    result = conv.WriteString(&mol,true);
  return result;
}

void linear_arrow_sort(vector<arrow_t> &arrows)
{
  point_t start;
  start.x=0;
  start.y=0;
  vector<arrow_t> new_arrows;
  while (!arrows.empty())
    {
      double d=FLT_MAX;
      int closest=0;
      bool found = false;
      for (int i=0; i<arrows.size(); i++)
	{
	  bool linebreak = false;
	  if (!new_arrows.empty())
	    linebreak = arrows[i].head.x>arrows[i].tail.x && new_arrows.back().head.x>new_arrows.back().tail.x 
	      && abs(arrows[i].head.y-arrows[i].tail.y)<5  && abs(new_arrows.back().head.y-new_arrows.back().tail.y)<5
	      && (min(arrows[i].head.y,arrows[i].tail.y)-max(new_arrows.back().head.y,new_arrows.back().tail.y)>MAX_FONT_HEIGHT ||
		  min(new_arrows.back().head.y,new_arrows.back().tail.y)-max(arrows[i].head.y,arrows[i].tail.y)>MAX_FONT_HEIGHT);
	      //&& arrows[i].tail.x<new_arrows.back().head.x;
	  
	  if (distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y)<d && (!linebreak || new_arrows.back().linebreak))
	    {
	      d = distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y);
	      closest = i;
	      found = true;
	    }
	}
     
      if (found)
	{
	  new_arrows.push_back(arrows[closest]);
	  start=arrows[closest].head;
	  arrows.erase(arrows.begin()+closest);
	}
      else
	{
	  start.x = 0;
	  start.y = 0;
	  if (!new_arrows.empty())
	    {
	      new_arrows.back().linebreak = true;
	      start.y = max(new_arrows.back().tail.y,new_arrows.back().head.y);
	    }
	}
    }
  arrows = new_arrows;
}

void check_the_last_arrow_linebreak(vector<arrow_t> &arrows,const vector<box_t> &page_of_boxes)
{
  if (arrows.back().head.x>arrows.back().tail.x && abs(arrows.back().head.y-arrows.back().tail.y)<5)
    {
      bool linebreak = true;
      for (int i=0; i<page_of_boxes.size(); i++)
	if (page_of_boxes[i].x1 > arrows.back().head.x && page_of_boxes[i].x1 - arrows.back().head.x < MAX_DISTANCE_BETWEEN_ARROWS
	    && page_of_boxes[i].y1 < arrows.back().head.y && page_of_boxes[i].y2 > arrows.back().head.y)
	  linebreak = false;
      arrows.back().linebreak = linebreak;
    }
}

void mark_reversible(vector<arrow_t> &arrows)
{
  if (arrows.size()<2) return;
  for (int i=0; i<arrows.size(); i++)
    for (int j=i+1; j<arrows.size(); j++)
      if (angle4(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,
		 arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y)<-D_T_TOLERANCE)
      {
	double d1=distance_from_bond_y(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,arrows[j].tail.x,arrows[j].tail.y);
	double d2=distance_from_bond_y(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,arrows[j].head.x,arrows[j].head.y);
	if (max(fabs(d1),fabs(d2))<2*MAX_FONT_HEIGHT)
	  {
	    
	    double l = distance(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y);
	    double l1 = distance_from_bond_x_a(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,arrows[j].tail.x,arrows[j].tail.y);
	    double l2 = distance_from_bond_x_a(arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y,arrows[j].head.x,arrows[j].head.y);
	    if (fabs(l-l1)<5 && fabs(l2)<5)
	      {
		double x1=0,y1=0, x2=1,y2=1;
		if (i>0)
		  {
		    x1=arrows[i-1].tail.x;
		    y1=arrows[i-1].tail.y;
		    int k = i+1;
		    while (k<arrows.size() && k==j) k++;
		    if (k==arrows.size())
		      {
			x2=arrows[i-1].head.x;
			y2=arrows[i-1].head.y;
		      }
		    else 
		      {
			x1=arrows[i-1].head.x;
			y1=arrows[i-1].head.y;
			x2=arrows[k].tail.x;
			y2=arrows[k].tail.y;
		      }
		  }
		double a1 = angle4(x1,y1,x2,y2,arrows[i].tail.x,arrows[i].tail.y,arrows[i].head.x,arrows[i].head.y);
		double a2 = angle4(x1,y1,x2,y2,arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y);
		if (a1>a2)
		  {
		    arrows[i].reversible = true;
		    arrows[j].remove=true;
		  }
		else
		  {
		    arrows[j].reversible = true;
		    arrows[i].remove=true;
		  }
	      }
	  }
      }
      
  vector<arrow_t>::iterator i=arrows.begin();
  while (i!=arrows.end())
    {
      if (i->remove) i=arrows.erase(i);
      else i++;
    }
}

double distance_from_box(const point_t &p, const box_t &b)
{
  if (p.x <= b.x1 && p.y <= b.y1)  return(distance(p.x,p.y,b.x1,b.y1));
  if (p.x > b.x1 && p.x < b.x2 && p.y < b.y1)  return(b.y1-p.y);
  if (p.x >= b.x2 && p.y <= b.y1)  return(distance(p.x,p.y,b.x2,b.y1));
  if (p.x < b.x1 && p.y > b.y1 && p.y < b.y2)  return(b.x1-p.x);
  if (p.x > b.x2 && p.y > b.y1 && p.y < b.y2)  return(p.x-b.x2);
if (p.x <= b.x1 && p.y >= b.y2)  return(distance(p.x,p.y,b.x1,b.y2));
if (p.x > b.x1 && p.x < b.x2 && p.y > b.y2)  return(p.y-b.y2);
if (p.x >= b.x2 && p.y >= b.y2)  return(distance(p.x,p.y,b.x2,b.y2));
return 0;
}

double distance_between_boxes(const box_t &a, const box_t &b)
{
  point_t p;
  p.x = a.x1;
  p.y = a.y1;
  double dab = distance_from_box(p,b);
  p.x = a.x1;
  p.y = a.y2;
 dab = min(dab,distance_from_box(p,b));
  p.x = a.x2;
  p.y = a.y1;
  dab = min(dab,distance_from_box(p,b));
  p.x = a.x2;
  p.y = a.y2;
  dab = min(dab,distance_from_box(p,b));

  p.x = b.x1;
  p.y = b.y1;
  dab = min(dab,distance_from_box(p,a));
  p.x = b.x2;
  p.y = b.y1;
  dab = min(dab,distance_from_box(p,a));
  p.x = b.x1;
  p.y = b.y2;
  dab = min(dab,distance_from_box(p,a));
  p.x = b.x2;
  p.y = b.y2;
  dab = min(dab,distance_from_box(p,a));
  return(dab);
}

void  sort_boxes_one_by_one(vector<int> &b, const vector<int> &a, const vector<box_t> &page_of_boxes, bool behind)
{
  box_t p=page_of_boxes[b[0]];

  set<int> existing;
  existing.insert(b[0]);
  if (behind)
    existing.insert(a[0]);
  else
    existing.insert(a.back());
  bool found = true;
  while (found)
    {
      found  = false;
      double d = FLT_MAX;
      int min_j = 0;    
      for (int j=0; j<page_of_boxes.size(); j++)
	if (d > distance_between_boxes(p,page_of_boxes[j]) && existing.find(j) == existing.end())
	  {
	    d = distance_between_boxes(p,page_of_boxes[j]);
	    min_j = j;
	  }
      if (d<MAX_DISTANCE_BETWEEN_ARROWS)
	{
	  b.push_back(min_j);
	  p = page_of_boxes[min_j];
	  found = true;
	  existing.insert(min_j);
	}
    }
}

void sort_boxes_from_arrows(vector < vector<int> > &before, const vector < vector<int> > &after,const vector<box_t> &page_of_boxes, bool behind)
{
 
  for (int i=0; i<before.size(); i++)
    if (!before[i].empty() && !after[i].empty())
    {
      sort_boxes_one_by_one(before[i],after[i],page_of_boxes,behind);
      if (behind)
	reverse(before[i].begin(),before[i].end());
    }
}

void arrange_structures_between_arrows_before(vector<arrow_t> &arrows,  vector < vector<int> > &before, const vector<box_t> &page_of_boxes, const vector<string> &page_of_structures)
{
  for (int j=0; j<arrows.size(); j++)
    {
      double rt = FLT_MAX;
      int i_min=0;

      for (int i=0; i<page_of_boxes.size(); i++)
	{
	  double r = distance_from_box(arrows[j].tail, page_of_boxes[i]);
	  double ry = distance_from_bond_y(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,(page_of_boxes[i].x2+page_of_boxes[i].x1)/2,(page_of_boxes[i].y2+page_of_boxes[i].y1)/2);
	  double l = distance(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y);
	  double rx1 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x1,page_of_boxes[i].y1);
	  double rx2 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x1,page_of_boxes[i].y2);
	  double rx3 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x2,page_of_boxes[i].y1);
	  double rx4 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x2,page_of_boxes[i].y2);
	  double cr = distance((page_of_boxes[i].x2+page_of_boxes[i].x1)/2,(page_of_boxes[i].y2+page_of_boxes[i].y1)/2,(arrows[j].tail.x+arrows[j].head.x)/2,(arrows[j].tail.y+arrows[j].head.y)/2);
	  if (rx1>=0 && rx1<=l && rx2>=0 && rx2<=l && rx3>=0 && rx3<=l && rx4>=0 && rx4<=l && 
	      cr<max(page_of_boxes[i].x2-page_of_boxes[i].x1, page_of_boxes[i].y2-page_of_boxes[i].y1))
	    {
	      string smi=convert_to_smiles_agent_structure(page_of_structures[i]);
		  if (!smi.empty())
			arrows[j].agent += " OSRA_AGENT_SMILES="+smi;
	    }
	  else if (fabs(ry)<min(page_of_boxes[i].x2-page_of_boxes[i].x1, page_of_boxes[i].y2-page_of_boxes[i].y1) && r<rt)
	    {
	      rt = r;
	      i_min = i;
	    }
	}

      if (rt<MAX_DISTANCE_BETWEEN_ARROWS)
	before[j].push_back(i_min);
    }
}

void arrange_structures_between_arrows_after(vector<arrow_t> &arrows,  vector < vector<int> > &after, const vector < vector<int> > &before,const vector<box_t> &page_of_boxes, const vector<string> &page_of_structures)
{
  for (int j=0; j<arrows.size(); j++)
    {
      double rh = FLT_MAX;
      int i_min=0;
      int previous_top = INT_MAX;
      int previous_bottom = INT_MAX;
      point_t p=arrows[j].head;
      if (arrows[j].linebreak)	p.x=0;
      if (before[j].empty()) continue;
      box_t box_before = page_of_boxes[before[j][0]];

      for (int i=0; i<page_of_boxes.size(); i++)
	{
	  double r = distance_from_box(p, page_of_boxes[i]);
	  double ry = distance_from_bond_y(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,(page_of_boxes[i].x2+page_of_boxes[i].x1)/2,(page_of_boxes[i].y2+page_of_boxes[i].y1)/2);
	  double l = distance(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y);
	  double rx1 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x1,page_of_boxes[i].y1);
	  double rx2 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x1,page_of_boxes[i].y2);
	  double rx3 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x2,page_of_boxes[i].y1);
	  double rx4 = distance_from_bond_x_a(arrows[j].tail.x,arrows[j].tail.y,arrows[j].head.x,arrows[j].head.y,page_of_boxes[i].x2,page_of_boxes[i].y2);
	  double cr = distance((page_of_boxes[i].x2+page_of_boxes[i].x1)/2,(page_of_boxes[i].y2+page_of_boxes[i].y1)/2,(arrows[j].tail.x+arrows[j].head.x)/2,(arrows[j].tail.y+arrows[j].head.y)/2);
	  bool agent_structure = rx1>=0 && rx1<=l && rx2>=0 && rx2<=l && rx3>=0 && rx3<=l && rx4>=0 && rx4<=l && cr<max(page_of_boxes[i].x2-page_of_boxes[i].x1, page_of_boxes[i].y2-page_of_boxes[i].y1);
	  if (!agent_structure)
	    {
	      if (arrows[j].linebreak)
		{
		  if (page_of_boxes[i].y1>box_before.y2 && page_of_boxes[i].y1-box_before.y2<MAX_DISTANCE_BETWEEN_ARROWS 
		  && ((r<rh && page_of_boxes[i].y1<previous_bottom) || page_of_boxes[i].y2<previous_top))
		    {
		      rh = r;
		      i_min = i;
		      previous_top = page_of_boxes[i].y1;
		      previous_bottom = page_of_boxes[i].y2;
		    }
		}
	      else if  (fabs(ry)<min(page_of_boxes[i].x2-page_of_boxes[i].x1, page_of_boxes[i].y2-page_of_boxes[i].y1)  && r<rh)
		{
		  rh = r;
		  i_min = i;
		}
	    }
	}

      if (arrows[j].linebreak)
	{
	  if (rh<FLT_MAX && previous_top<INT_MAX && page_of_boxes[i_min].y1>box_before.y2 && page_of_boxes[i_min].y1-box_before.y2<MAX_DISTANCE_BETWEEN_ARROWS)
	    after[j].push_back(i_min);
	}
      else if (rh<MAX_DISTANCE_BETWEEN_ARROWS)
	after[j].push_back(i_min);
    }
}

void arrange_plus_sings_between_boxes(const vector < vector<int> > &before,const vector<box_t> &page_of_boxes,const vector<plus_t> &pluses, vector < vector <int> > &is_plus)
{
  for (int i=0; i<before.size(); i++)
    for (int j=1; j<before[i].size(); j++)
      {
	int l = before[i][j-1];
	int k = before[i][j];
	if (k>=0 && l>=0)
	  {
	    box_t a = page_of_boxes[l];
	    box_t b = page_of_boxes[k];
	    for (int m=0; m<pluses.size(); m++)
	      {
		//double d=distance_from_bond_y((a.x1+a.x2)/2,(a.y1+a.y2)/2,(b.x1+b.x2)/2,(b.y1+b.y2)/2,pluses[m].x, pluses[m].y);
		if (pluses[m].center.y>max(a.y1,b.y1) && pluses[m].center.y<min(a.y2,b.y2) && ((pluses[m].center.x>a.x2 && pluses[m].center.x<b.x1) ||  (pluses[m].center.x>b.x2 && pluses[m].center.x<a.x1)))
		  {
		    is_plus[k][l] = m;
		    is_plus[l][k] = m;
		  }
		// after plus things can be on the next line
		//d = pluses[m].y - (a.y2 + a.y1)/2;
		if (pluses[m].center.x>a.x2 &&  pluses[m].center.y>a.y1 && pluses[m].center.y<a.y2 && b.y1>a.y2 && pluses[m].center.x-a.x2<MAX_DISTANCE_BETWEEN_ARROWS && b.y1-a.y2<MAX_DISTANCE_BETWEEN_ARROWS)
		  {
		    is_plus[k][l] = m;
		    is_plus[l][k] = m;
		  }
		/*if (pluses[m].x>b.x2 &&  pluses[m].y>b.y1 && pluses[m].y<b.y2 && a.y1>b.y2 && pluses[m].x-b.x2<MAX_DISTANCE_BETWEEN_ARROWS && a.y1-b.y2<MAX_DISTANCE_BETWEEN_ARROWS)
		  {
		    is_plus[k][l] = true;
		    is_plus[l][k] = true;
		  }*/
		
	      }
	  }
      }
}


void arrange_reactions(vector<arrow_t> &arrows, const vector<box_t> &page_of_boxes, const vector<plus_t> &pluses, vector<string> &results, vector<box_t> &rbox,
		       const vector<string> &page_of_structures,  const string &output_format)
{
  vector < vector<int> > before;
  vector < vector<int> > after;
  if (arrows.empty() || page_of_boxes.empty()) return;
  // Find average distance between nearest arrows and standard deviation
  vector<arrow_t> arrows_by_closest;
  vector<double> dist_arrows;
  double avg_dist_arrow=0, avg_dist_arrow_squared=0, std_dev_arrow=0;
  int n_arrows = arrows.size();
  point_t start;
  start.x=0;
  start.y=0;

  while (!arrows.empty())
    {
      double d = FLT_MAX;
      int min_i=0;
      for (int i=0; i<arrows.size(); i++)
	if (d>min(distance(start.x,start.y,arrows[i].head.x,arrows[i].head.y),distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y)))
	  {
	    d = min(distance(start.x,start.y,arrows[i].head.x,arrows[i].head.y),distance(start.x,start.y,arrows[i].tail.x,arrows[i].tail.y));
	    min_i = i;
	  }
      arrows_by_closest.push_back(arrows[min_i]);
      dist_arrows.push_back(d);
      start = arrows[min_i].head;
      avg_dist_arrow += d;
      avg_dist_arrow_squared +=d*d;
      arrows.erase(arrows.begin()+min_i);
    }
  avg_dist_arrow /= n_arrows;
  avg_dist_arrow_squared /= n_arrows;
  std_dev_arrow = sqrt(avg_dist_arrow_squared - avg_dist_arrow*avg_dist_arrow);

  mark_reversible(arrows_by_closest);

  // Break arrows into close-knit groups
  vector < vector <arrow_t> > arrow_groups;
  int i=0;
  while (i<n_arrows)
    {
      vector <arrow_t> group;
      group.push_back(arrows_by_closest[i]);
      i++;
      while (i<n_arrows && dist_arrows[i]<avg_dist_arrow+2*std_dev_arrow)
	{
	  group.push_back(arrows_by_closest[i]);
	  i++;
	}
      arrow_groups.push_back(group);
    }

  // arrange arrows in head to tail fashion
  for (int i=0; i<arrow_groups.size(); i++)
    {
      linear_arrow_sort(arrow_groups[i]);
      check_the_last_arrow_linebreak(arrow_groups[i],page_of_boxes);
    }
  // combine groups into a flat list of arrows
  arrows.clear();
  for (int i=0; i<arrow_groups.size(); i++)
    {
      // arrange structures to best fit between arrows
      vector < vector<int> > before_group(arrow_groups[i].size());
      arrange_structures_between_arrows_before(arrow_groups[i],before_group,page_of_boxes,page_of_structures);
      vector < vector<int> > after_group(arrow_groups[i].size());
      arrange_structures_between_arrows_after(arrow_groups[i],after_group,before_group,page_of_boxes,page_of_structures);
    


      sort_boxes_from_arrows(before_group,after_group,page_of_boxes,true);
      sort_boxes_from_arrows(after_group,before_group,page_of_boxes,false);
      /*for (int ii=0;ii<before_group.size(); ii++)
   {
     for (int j=0; j<before_group[ii].size(); j++)
    	cout<<before_group[ii][j]<<" ";
     cout<<">>> ";
     for (int j=0; j<after_group[ii].size(); j++)
    	cout<<after_group[ii][j]<<" ";
     cout<<endl;
   }
      */
      for (int j=0; j<arrow_groups[i].size(); j++)
	arrows.push_back(arrow_groups[i][j]);

      for (int j=0; j<before_group.size(); j++)
	before.push_back(before_group[j]);

      for (int j=0; j<after_group.size(); j++)
	after.push_back(after_group[j]);
    }

  //for (int i=0; i<arrows.size(); i++)
    //cout<<arrows[i].tail.x<<","<<arrows[i].tail.y<<" "<<arrows[i].head.x<<","<<arrows[i].head.y<<" "<<arrows[i].linebreak<<endl;
  //for (int i=0; i<pluses.size(); i++)
    //    cout<<pluses[i].x<<" "<<pluses[i].y<<endl;
  /*for (int ii=0;ii<before.size(); ii++)
   {
     for (int j=0; j<before[ii].size(); j++)
    	cout<<before[ii][j]<<" ";
     cout<<">>> ";
     for (int j=0; j<after[ii].size(); j++)
    	cout<<after[ii][j]<<" ";
     cout<<endl;
     }*/
  
 
  vector < vector <int> > is_plus_before(page_of_boxes.size(), vector <int> (page_of_boxes.size(), -1));
  arrange_plus_sings_between_boxes(before,page_of_boxes,pluses, is_plus_before);
 vector < vector <int> > is_plus_after(page_of_boxes.size(), vector <int> (page_of_boxes.size(), -1));
  arrange_plus_sings_between_boxes(after,page_of_boxes,pluses, is_plus_after);

/*   for (int ii=0;ii<before.size(); ii++)
   {
     for (int j=0; j<before[ii].size(); j++)
       {
	 cout<<before[ii][j];
	 if (j<before[ii].size()-1 && is_plus_before[before[ii][j]][before[ii][j+1]])
	   cout<<"+";
	 else
	   cout<<" ";
       }
     cout<<">>> ";
     for (int j=0; j<after[ii].size(); j++)
       {
	 cout<<after[ii][j];
	 if (j<after[ii].size()-1 && is_plus_after[after[ii][j]][after[ii][j+1]])
	   cout<<"+";
	 else
	   cout<<" ";
       }

     cout<<endl;
     }*/

  // extract reactions, if any
  for (int i=0; i<arrows.size(); i++)
    {
      vector <int> r,p;
      box_t box;
      box.x1 = arrows[i].min_x;
      box.y1 = arrows[i].min_y;
      box.x2 = arrows[i].max_x;
      box.y2 = arrows[i].max_y;

      int ii = before[i].size()-1;
      while (ii>=0 && before[i][ii]<0) ii--;

      if (ii>=0 && before[i][ii]>=0)
	{
	  int k = before[i][ii];
	  r.push_back(k);
	  box.x1 = min(box.x1,page_of_boxes[k].x1);
	  box.y1 = min(box.y1,page_of_boxes[k].y1);
	  box.x2 = max(box.x2,page_of_boxes[k].x2);
	  box.y2 = max(box.y2,page_of_boxes[k].y2);
	}

      for (int j=ii-1; j>=0; j--)
	{
	  int k = before[i][j];
	  int l = before[i][j+1];
	  if (k>=0 && l>=0)
	    {
	      int m = is_plus_before[k][l];
	      if (m>=0)
		{
		  r.push_back(k);
		  box.x1 = min(box.x1,page_of_boxes[k].x1);
		  box.y1 = min(box.y1,page_of_boxes[k].y1);
		  box.x2 = max(box.x2,page_of_boxes[k].x2);
		  box.y2 = max(box.y2,page_of_boxes[k].y2);

		  box.x1 = min(box.x1,pluses[m].min_x);
		  box.y1 = min(box.y1,pluses[m].min_y);
		  box.x2 = max(box.x2,pluses[m].max_x);
		  box.y2 = max(box.y2,pluses[m].max_y);
		}
	      else
		break;
	    }
	}

      if (!after[i].empty())
	{
	  ii = 0;
	  while (ii<after[i].size() && after[i][ii]<0) ii++;

	  if (ii<after[i].size() && after[i][ii]>=0)
	    {
	      int k = after[i][ii];
	      p.push_back(k);
	      box.x1 = min(box.x1,page_of_boxes[k].x1);
	      box.y1 = min(box.y1,page_of_boxes[k].y1);
	      box.x2 = max(box.x2,page_of_boxes[k].x2);
	      box.y2 = max(box.y2,page_of_boxes[k].y2);
	    }
	  for (int j=ii+1; j<after[i].size(); j++)
	    {
	      int k = after[i][j];
	      int l = after[i][j-1];
	      if (k>=0 && l>=0)
		{
		  int m = is_plus_after[k][l];
		  if (m>=0)
		    {
		      p.push_back(k);
		      box.x1 = min(box.x1,page_of_boxes[k].x1);
		      box.y1 = min(box.y1,page_of_boxes[k].y1);
		      box.x2 = max(box.x2,page_of_boxes[k].x2);
		      box.y2 = max(box.y2,page_of_boxes[k].y2);

		      box.x1 = min(box.x1,pluses[m].min_x);
		      box.y1 = min(box.y1,pluses[m].min_y);
		      box.x2 = max(box.x2,pluses[m].max_x);
		      box.y2 = max(box.y2,pluses[m].max_y);
		    }
		  else
		    break;
		}
	    }
	}

      if (!r.empty() && !p.empty())
	{
	  string result=convert_page_to_reaction(page_of_structures,output_format, r, p, arrows[i].agent,arrows[i].reversible);
	  trim(result);
	  if (!result.empty())
	    {
	      results.push_back(result);
	      rbox.push_back(box);
	    }
	}
    }
}
