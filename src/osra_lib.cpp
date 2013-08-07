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

#include <stdio.h> // fclose
#include <stdlib.h> // malloc(), free()
#include <math.h> // fabs(double)
#include <float.h> // FLT_MAX
#include <limits.h> // INT_MAX

#include <list> // sdt::list
#include <vector> // std::vector
#include <set>
#include <algorithm> // std::sort, std::min(double, double), std::max(double, double)
#include <iostream> // std::ostream, std::cout
#include <fstream> // std::ofstream, std::ifstream
#include <sstream> // std:ostringstream

#include <Magick++.h>

extern "C" {
#include <potracelib.h>
#include <pgm2asc.h>
}

#include <openbabel/oberror.h>

#include "osra.h"
#include "osra_grayscale.h"
#include "osra_segment.h"
#include "osra_fragments.h"
#include "osra_labels.h"
#include "osra_thin.h"
#include "osra_common.h"
#include "osra_structure.h"
#include "osra_lib.h"
#include "osra_ocr.h"
#include "osra_openbabel.h"
#include "osra_reaction.h"
#include "osra_anisotropic.h"
#include "osra_stl.h"
#include "unpaper.h"
#include "config.h" // DATA_DIR

using namespace std;
using namespace Magick;

void set_select_resolution(vector<int>  &select_resolution, int input_resolution)
{
  if (input_resolution == 0)
    {
      select_resolution[0] = 72;
      select_resolution[1] = 150;
      select_resolution[2] = 300;
      select_resolution[3] = 300;  // No thinning
      select_resolution[4] = 500;
    }
}

double set_threshold(double threshold,int resolution)
{
  double THRESHOLD_BOND = threshold;

  if (THRESHOLD_BOND < 0.0001)
    {
      if (resolution >= 150)
        {
          THRESHOLD_BOND = THRESHOLD_GLOBAL;
        }
      else
        {
          THRESHOLD_BOND = THRESHOLD_LOW_RES;
        }
    }
  return THRESHOLD_BOND;
}

int load_superatom_spelling_maps(map<string, string> &spelling,map<string, string> &superatom, const string &osra_dir,
                                 const string &spelling_file, const string &superatom_file, bool verbose)
{
// Loading the program data files into maps:
  if (!((spelling_file.length() != 0 && load_config_map(spelling_file, spelling))
        || load_config_map(string(DATA_DIR) + "/" + SPELLING_TXT, spelling) || load_config_map(osra_dir + "/" + SPELLING_TXT, spelling)))
    {
      cerr << "Cannot open " << SPELLING_TXT << " file (tried locations \"" << DATA_DIR << "\", \"" << osra_dir
           << "\"). Specify the custom file location via -l option." << endl;
      return ERROR_SPELLING_FILE_IS_MISSING;
    }

  if (!((superatom_file.length() != 0 && load_config_map(superatom_file, superatom))
        || load_config_map(string(DATA_DIR) + "/" + SUPERATOM_TXT, superatom) || load_config_map(osra_dir + "/"
            + SUPERATOM_TXT, superatom)))
    {
      cerr << "Cannot open " << SUPERATOM_TXT << " file (tried locations \"" << DATA_DIR << "\", \"" << osra_dir
           << "\"). Specify the custom file location via -a option." << endl;
      return ERROR_SUPERATOM_FILE_IS_MISSING;
    }

  if (verbose)
    cout << "spelling (size: " << spelling.size() << ") and superatom (size: " << superatom.size() << ") dictionaries are loaded." << endl;
  return 0;
}

void create_thick_box(Image &orig_box,Image &thick_box,int &width,int &height,int &resolution,int &working_resolution,double &box_scale,
                      ColorGray bgColor, double THRESHOLD_BOND, int res_iter, bool &thick, bool jaggy)
{
  if (resolution >= 300)
    {
      int max_hist;
      double nf45;
      double nf =
        noise_factor(orig_box, width, height, bgColor, THRESHOLD_BOND, resolution, max_hist, nf45);

      //if (max_hist < 5) thick = false;
      if (res_iter == NUM_RESOLUTIONS-2)  // no thinning
	{
	  working_resolution = 300;
	  thick_box = orig_box;
	  width = thick_box.columns();
	  height = thick_box.rows();
	  thick = false;
	  //nf = noise_factor(orig_box, width, height, bgColor, THRESHOLD_BOND, resolution, max_hist, nf45);
	}
      if (res_iter == NUM_RESOLUTIONS-1)
        {
          if (max_hist >= 6)
            {
              int new_resolution = max_hist * 300 / 4;
              int percent = (100 * 300) / new_resolution;
              resolution = new_resolution;
              ostringstream scale;
              scale << percent << "%";
              orig_box.scale(scale.str());
              box_scale /= (double) percent/100;
              working_resolution = 300;
              thick_box = orig_box;
              width = thick_box.columns();
              height = thick_box.rows();
	      nf = noise_factor(orig_box, width, height, bgColor, THRESHOLD_BOND, resolution, max_hist, nf45);
            }
          else
            {
              resolution = 500;
              int percent = (100 * 300) / resolution;
              ostringstream scale;
              scale << percent << "%";
              orig_box.scale(scale.str());
              box_scale /= (double) percent/100;
              working_resolution = 300;
              thick_box = orig_box;
              width = thick_box.columns();
              height = thick_box.rows();
              thick = false;
              nf = noise_factor(orig_box, width, height, bgColor, THRESHOLD_BOND, resolution,  max_hist, nf45);
            }
        }
      if (jaggy)
        {
          orig_box.scale("50%");
          box_scale *= 2;
          working_resolution = 150;
          //orig_box.scale("33%");
          //box_scale *= 3;
          // working_resolution = 100;
          thick_box = orig_box;
          width = thick_box.columns();
          height = thick_box.rows();
        }
      else if (nf > 0.5 && nf < 1. && max_hist <= 6)// && res_iter != 3 && max_hist <= 6)
        try
          {
            thick_box = anisotropic_smoothing(orig_box, width, height, 20, 0.3, 1.0, 0.6, 2);
          }
        catch (...)
          {
            thick_box = orig_box;
          }
      /*else if (nf45 > 0.9 && nf45 < 1.2 && max_hist == 3)
      {
      //orig_box = anisotropic_smoothing(thick_box, width, height, 60, 0.3, 0.6, 4., 2.);
      orig_box.scale("50%");
      thick_box = orig_box;
      //working_resolution = 150;
      width = thick_box.columns();
      height = thick_box.rows();
      //thick = false;
      }*/
      else
        thick_box = orig_box;
    }
  else if (resolution < 300 && resolution > 150)
    {
      int nw = width * 300 / resolution;
      int nh = height * 300 / resolution;
      thick_box = anisotropic_scaling(orig_box, width, height, nw, nh);
      width = thick_box.columns();
      height = thick_box.rows();
      int percent = (100 * 300) / resolution;
      ostringstream scale;
      scale << percent << "%";
      orig_box.scale(scale.str());
      box_scale /= (double) percent/100;
      working_resolution = 300;
    }
  else
    thick_box = orig_box;
}

potrace_state_t * const  raster_to_vector(Image &box,ColorGray bgColor, double THRESHOLD_BOND,int width,int height,int working_resolution)
{
  potrace_param_t * const param = potrace_param_default();
  param->alphamax = 5e-324; // this has been changed in potrace-1.11 
  //param->turnpolicy = POTRACE_TURNPOLICY_MINORITY;
  param->turdsize = 0;

  param->turnpolicy = POTRACE_TURNPOLICY_MINORITY;
  double c_width = 1. * width * 72 / working_resolution;
  double c_height = 1. * height * 72 / working_resolution;
  if (c_height * c_width < SMALL_PICTURE_AREA)
    param->turnpolicy = POTRACE_TURNPOLICY_BLACK;

  potrace_bitmap_t * const bm = bm_new(width, height);
  for (int i = 0; i < width; i++)
    for (int j = 0; j < height; j++)
      BM_PUT(bm, i, j, get_pixel(box, bgColor, i, j, THRESHOLD_BOND));

  potrace_state_t * const st = potrace_trace(param, bm);
  if (bm != NULL)
    {
      free(bm->map);
      free(bm);
    }
  potrace_param_free(param);
  return(st);
}

void rotate_point(int &x, int &y, int midX, int midY, double rotation)
{
// create 2D rotation matrix
  float sinval = sin(rotation); // no use of sincos()-function for compatibility, no performace bottleneck anymore anyway
  float cosval = cos(rotation);
  float m11 = cosval;
  float m12 = sinval;
  float m21 = -sinval;
  float m22 = cosval;

  int dX = x - midX;
  int dY = y - midY;

  int diffX = dX * m11 + dY * m21;
  int diffY = dX * m12 + dY * m22;

  x = midX+diffX;
  y = midY+diffY;
}

void rotate_coordinate_box(box_t &coordinate_box,double rotation,int width,int height)
{
  int midX = width/2;
  int midY = height/2;
  int x1 = coordinate_box.x1;
  int y1 = coordinate_box.y1;
  int x2 = coordinate_box.x2;
  int y2 = coordinate_box.y2;

  rotate_point(x1,y1,midX,midY,rotation);
  rotate_point(x2,y1,midX,midY,rotation);
  rotate_point(x1,y2,midX,midY,rotation);
  rotate_point(x2,y2,midX,midY,rotation);

  coordinate_box.x1 = min(x1,x2);
  coordinate_box.x2 = max(x1,x2);
  coordinate_box.y1 = min(y1,y2);
  coordinate_box.y2 = max(y1,y2);
}

void split_fragments_and_assemble_structure_record(vector<atom_t> &atom,
						   int n_atom, 
						   vector<bond_t>  &bond, 
						   int n_bond, 
						   const vector<box_t> &boxes,
						   int l,int k,
						   int resolution,
						   int res_iter, 
						   const string &output_image_file_prefix,
						   Image &image,
						   Image &orig_box,
						   int real_font_width,int real_font_height,
						   double thickness, 
						   double avg_bond_length,
						   const map<string, string> &superatom,
						   int real_atoms, int real_bonds, 
						   int bond_max_type,
						   double box_scale, double page_scale, double rotation, int unpaper_dx, int unpaper_dy,
						   string output_format,
						   const string &embedded_format,
						   bool is_reaction,
						   bool show_confidence,
						   bool show_resolution_guess,
						   bool show_page,
						   bool show_coordinates,
						   bool show_avg_bond_length,
						   vector<vector<string> > &array_of_structures,
						   vector<vector<double> > &array_of_avg_bonds,
						   vector<vector<double> > &array_of_ind_conf,
						   vector<vector<Image> > &array_of_images,
						   vector<vector<box_t> > &array_of_boxes,
						   int &total_boxes,
						   double &total_confidence,
						   int n_letters,
						   bool show_learning,
						   int resolution_iteration,
						   bool verbose)
{
  vector<atom_t> frag_atom;
  vector<bond_t> frag_bond;

  if (real_atoms > MIN_A_COUNT && real_atoms < MAX_A_COUNT && real_bonds < MAX_A_COUNT && bond_max_type>0 && bond_max_type<5)
    {
      int num_frag;
      num_frag = resolve_bridge_bonds(atom, n_atom, bond, n_bond, 2 * thickness, avg_bond_length, superatom, verbose);
      collapse_bonds(atom, bond, n_bond, avg_bond_length / 4);
      collapse_atoms(atom, bond, n_atom, n_bond, 3);
      remove_zero_bonds(bond, n_bond, atom);
      extend_terminal_bond_to_bonds(atom, bond, n_bond, avg_bond_length, 7, 0);


      remove_small_terminal_bonds(bond, n_bond, atom, avg_bond_length);
      n_bond = reconnect_fragments(bond, n_bond, atom, avg_bond_length);

      collapse_atoms(atom, bond, n_atom, n_bond, 1);
      mark_terminal_atoms(bond, n_bond, atom, n_atom);
      const vector<vector<int> > &frags = find_fragments(bond, n_bond, atom);
      vector<fragment_t> fragments = populate_fragments(frags, atom);
      std::sort(fragments.begin(), fragments.end(), comp_fragments);
      for (unsigned int i = 0; i < fragments.size(); i++)
        {
          if (verbose)
            cout << "Considering fragment #" << i + 1 << " " << fragments[i].x1 << "x" << fragments[i].y1 << "-" << fragments[i].x2 << "x"
                 << fragments[i].y2 << ", atoms: " << fragments[i].atom.size() << '.' << endl;

          if (fragments[i].atom.size() > MIN_A_COUNT)
            {
              frag_atom.clear();
              for (int a = 0; a < n_atom; a++)
                {
                  frag_atom.push_back(atom[a]);
                  frag_atom[a].exists = false;
                }

              for (unsigned int j = 0; j < fragments[i].atom.size(); j++)
                frag_atom[fragments[i].atom[j]].exists = atom[fragments[i].atom[j]].exists;

              frag_bond.clear();
              for (int b = 0; b < n_bond; b++)
                {
                  frag_bond.push_back(bond[b]);
                }

              remove_zero_bonds(frag_bond, n_bond, frag_atom);

              double confidence = 0;
              molecule_statistics_t molecule_statistics;
              int page_number = l + 1;
              box_t coordinate_box,rel_box;
	       if (fragments.size()>1)
		{
		  coordinate_box.x1 = (int) (-(double)page_scale * unpaper_dx + (double) page_scale * boxes[k].x1 + (double) page_scale * box_scale * fragments[i].x1 - (double) page_scale * FRAME);
		  coordinate_box.y1 = (int) (-(double)page_scale * unpaper_dy + (double) page_scale * boxes[k].y1 + (double) page_scale * box_scale * fragments[i].y1 - (double) page_scale * FRAME);
		  coordinate_box.x2 = (int) (-(double)page_scale * unpaper_dx + (double) page_scale * boxes[k].x1 + (double) page_scale * box_scale * fragments[i].x2 - (double) page_scale * FRAME);
		  coordinate_box.y2 = (int) (-(double)page_scale * unpaper_dy + (double) page_scale * boxes[k].y1 + (double) page_scale * box_scale * fragments[i].y2 - (double) page_scale * FRAME);
		  //rotate_coordinate_box(coordinate_box,rotation,image.columns(),image.rows());
		  rel_box.x1 = (int)((double)boxes[k].x1 + (double) box_scale * fragments[i].x1 - FRAME);
		  rel_box.y1 = (int)((double)boxes[k].y1 + (double) box_scale * fragments[i].y1 - FRAME);
		  rel_box.x2 = (int)((double)boxes[k].x1 + (double) box_scale * fragments[i].x2 - FRAME);
		  rel_box.y2 = (int)((double)boxes[k].y1 + (double) box_scale * fragments[i].y2 - FRAME);
		}
		else
		{
		  coordinate_box.x1 = (int) (-(double)page_scale * unpaper_dx + (double) page_scale * boxes[k].x1);
		  coordinate_box.y1 = (int) (-(double)page_scale * unpaper_dy + (double) page_scale * boxes[k].y1);
		  coordinate_box.x2 = (int) (-(double)page_scale * unpaper_dx + (double) page_scale * boxes[k].x2);
		  coordinate_box.y2 = (int) (-(double)page_scale * unpaper_dy + (double) page_scale * boxes[k].y2);
		  //rotate_coordinate_box(coordinate_box,rotation,image.columns(),image.rows());
		  rel_box.x1 = boxes[k].x1;
		  rel_box.y1 = boxes[k].y1;
		  rel_box.x2 = boxes[k].x2;
		  rel_box.y2 = boxes[k].y2;
		  }

              if (verbose)
                cout << "Coordinate box: " << coordinate_box.x1 << "x" << coordinate_box.y1 << "-" << coordinate_box.x2 << "x"
                     << coordinate_box.y2 << "." << endl;
	      if (is_reaction)
		output_format = SUBSTITUTE_REACTION_FORMAT;

              string structure =
                get_formatted_structure(frag_atom, frag_bond, n_bond, output_format, embedded_format,
                                        molecule_statistics, confidence,
                                        show_confidence, avg_bond_length, page_scale * box_scale * avg_bond_length,
                                        show_avg_bond_length,
                                        show_resolution_guess ? &resolution : NULL,
                                        show_page ? &page_number : NULL,
                                        show_coordinates ? &coordinate_box : NULL, superatom, n_letters, show_learning, resolution_iteration, verbose);

              if (molecule_statistics.fragments > 0 && molecule_statistics.fragments < MAX_FRAGMENTS && molecule_statistics.num_atoms>MIN_A_COUNT && molecule_statistics.num_bonds>0)
                {
                  array_of_structures[res_iter].push_back(structure);
                  array_of_avg_bonds[res_iter].push_back(page_scale * box_scale * avg_bond_length);
                  array_of_ind_conf[res_iter].push_back(confidence);
		  array_of_boxes[res_iter].push_back(rel_box);
                  total_boxes++;
                  total_confidence += confidence;
		  if (verbose)
		    cout<<"Result: "<<res_iter<<" "<<structure<<" "<<confidence<<endl;
                  if (!output_image_file_prefix.empty())
                    {
                      Image tmp = image;
		      if (!is_reaction)
			{
			  Geometry geometry =
			    (fragments.size() > 1) ? Geometry(box_scale * fragments[i].x2 - box_scale * fragments[i].x1, //
							      box_scale * fragments[i].y2 - box_scale * fragments[i].y1, //
							      boxes[k].x1 + box_scale * fragments[i].x1 - FRAME , //
							      boxes[k].y1 + box_scale * fragments[i].y1 - FRAME )
			    : Geometry(boxes[k].x2 - boxes[k].x1, boxes[k].y2 - boxes[k].y1, boxes[k].x1, boxes[k].y1);
			  
			  try
			    {
			      tmp.crop(geometry);
			    }
			  catch (...)
			    {
			      tmp = orig_box;
			    }
			}

                      array_of_images[res_iter].push_back(tmp);
                    }
                }
            }
        }
    }
}

int count_recognized_chars(vector<atom_t>  &atom, vector<bond_t>& bond)
{
  string char_filter = "oOcCNHsSBMeEXYZRPp23456789AF";
  set<int> atoms;
  for (int i=0; i<bond.size(); i++)
    if (bond[i].exists)
      {
	atoms.insert(bond[i].a);
	atoms.insert(bond[i].b);
      }
  int r = 0;
  for (set<int>::iterator a=atoms.begin(); a != atoms.end(); a++)
    for (int i=0; i<atom[*a].label.size(); i++)
      {
	if (char_filter.find(atom[*a].label[i]) != string::npos)
	  r++;
      }
  return r;
}

extern job_t *OCR_JOB;
extern job_t *JOB;

//int global_init_state;

// Function: osra_init()
//
// Initialises OSRA library. Should be called at e.g. program startup. This function is automatically called for both SO library and CLI utility.
// See this section for details about library init/cleanup: http://www.faqs.org/docs/Linux-HOWTO/Program-Library-HOWTO.html#INIT-AND-CLEANUP
// Below attribute marker is GNU compiler specific.
void __attribute__ ((constructor)) osra_init()
{
  // Necessary for GraphicsMagick-1.3.8 according to http://www.graphicsmagick.org/1.3/NEWS.html#january-21-2010:
  MagickLib::InitializeMagick(NULL);

  osra_ocr_init();

// global_init_state = osra_openbabel_init();

  //if (global_init_state != 0)
  //cerr << "OpenBabel initialization failure." << endl;

  srand(1);
}

// Function: osra_destroy()
//
// Releases all resources allocated by OSRA library. Should be called at e.g. program exit. This function is automatically called for both SO library and CLI utility.
// See this section for details about library init/cleanup: http://www.faqs.org/docs/Linux-HOWTO/Program-Library-HOWTO.html#INIT-AND-CLEANUP
// Below attribute marker is GNU compiler specific.
void __attribute__ ((destructor)) osra_destroy()
{
#ifdef OSRA_LIB
  MagickLib::DestroyMagick();
#endif
  osra_ocr_destroy();
}

int osra_process_image(
#ifdef OSRA_LIB
  const char *image_data,
  int image_length,
  ostream &structure_output_stream,
#else
  const string &input_file,
  const string &output_file,
#endif
  int rotate,
  bool invert,
  int input_resolution,
  double threshold,
  int do_unpaper,
  bool jaggy,
  bool adaptive_option,
  string output_format,
  string embedded_format,
  bool show_confidence,
  bool show_resolution_guess,
  bool show_page,
  bool show_coordinates,
  bool show_avg_bond_length,
  bool show_learning,
  const string &osra_dir,
  const string &spelling_file,
  const string &superatom_file,
  bool debug,
  bool verbose,
  const string &output_image_file_prefix,
  const string &resize
)
{
  //if (global_init_state != 0)
  // return global_init_state;

  std::transform(output_format.begin(), output_format.end(), output_format.begin(), ::tolower);
  std::transform(embedded_format.begin(), embedded_format.end(), embedded_format.begin(), ::tolower);

  map<string, string> spelling, superatom;
  int err = load_superatom_spelling_maps(spelling, superatom, osra_dir, spelling_file, superatom_file, verbose);
  if (err != 0) return err;

  string type;

#ifdef OSRA_LIB
  Blob blob(image_data, image_length);
#endif

  try
    {
      Image image_typer;
#ifdef OSRA_LIB
      image_typer.ping(blob);
#else
      image_typer.ping(input_file);
#endif
      type = image_typer.magick();
    }
  catch (...)
    {
      // Unfortunately, GraphicsMagick does not throw exceptions in all cases, so it behaves inconsistent, see
      // https://sourceforge.net/tracker/?func=detail&aid=3022955&group_id=40728&atid=428740
    }

  if (type.empty())
    {
#ifdef OSRA_LIB
      cerr << "Cannot detect blob image type" << endl;
#else
      cerr << "Cannot open file \"" << input_file << '"' << endl;
#endif
      return ERROR_UNKNOWN_IMAGE_TYPE;
    }

  if (verbose)
    cout << "Image type: " << type << '.' << endl;

#ifndef OSRA_LIB
  ofstream outfile;

  if (!output_file.empty())
    {
      outfile.open(output_file.c_str(), ios::out | ios::trunc);
      if (outfile.bad() || !outfile.is_open())
        {
          cerr << "Cannot open file \"" << output_file << "\" for output" << endl;
          return ERROR_OUTPUT_FILE_OPEN_FAILED;
        }
    }
#endif

  if (input_resolution == 0 && (type == "PDF" || type == "PS"))
    input_resolution = 300;

  if (show_coordinates && rotate != 0)
    {
      cerr << "Showing the box coordinates is currently not supported together with image rotation and is therefore disabled." << endl;
#ifdef OSRA_LIB
      return ERROR_ILLEGAL_ARGUMENT_COMBINATION;
#else
      show_coordinates = false;
#endif
    }

  if (!embedded_format.empty() && !(output_format == "sdf" && (embedded_format == "inchi" || embedded_format == "smi"
                                    || embedded_format == "can")))
    {
      cerr << "Embedded format option is only possible if output format is SDF and option can have only inchi, smi, or can values." << endl;
      return ERROR_ILLEGAL_ARGUMENT_COMBINATION;
    }
  bool is_reaction = false;
  if (output_format == "cmlr" || output_format == "rsmi" || output_format =="rxn")
    is_reaction = true;

#ifdef OSRA_LIB
  int page = 1;
#else
  int page = count_pages(input_file);
#endif

  vector<vector<string> > pages_of_structures(page, vector<string> (0));
  vector<vector<Image> > pages_of_images(page, vector<Image> (0));
  vector<vector<double> > pages_of_avg_bonds(page, vector<double> (0));
  vector<vector<double> > pages_of_ind_conf(page, vector<double> (0));
  vector<vector<box_t> > pages_of_boxes(page, vector<box_t> (0));
  vector<vector<arrow_t> > arrows(page, vector<arrow_t>(0));
  vector<vector<plus_t> > pluses(page, vector<plus_t>(0));
 

  int total_structure_count = 0;
  int num_resolutions = NUM_RESOLUTIONS;
  if (input_resolution != 0)
    num_resolutions = 1;
  vector<double> array_of_confidence(num_resolutions, 0);
  vector<int> boxes_per_res(num_resolutions,0);
  vector<int> select_resolution(num_resolutions, input_resolution);
  set_select_resolution(select_resolution,input_resolution);
  vector<vector<vector<string> > > array_of_structures_page(page,vector<vector<string> >(num_resolutions));
  vector<vector<vector<double> > > array_of_avg_bonds_page(page,vector<vector<double> >(num_resolutions)), array_of_ind_conf_page(page,vector<vector<double> >(num_resolutions));
  vector<vector<vector<Image> > > array_of_images_page(page,vector<vector<Image> > (num_resolutions));
  vector<vector<vector<box_t> > > array_of_boxes_page(page,vector<vector<box_t> >(num_resolutions));

  #pragma omp parallel for default(shared) private(OCR_JOB,JOB)
  for (int l = 0; l < page; l++)
    {
      Image image;
      double page_scale=1;
      
      int ttt = 0;

      if (verbose)
        cout << "Processing page " << (l+1) << " out of " << page << "..." << endl;

      ostringstream density;
      density << input_resolution << "x" << input_resolution;
      image.density(density.str());

      if (type == "PDF" || type == "PS")
        page_scale *= (double) 72 / input_resolution;

#ifdef OSRA_LIB
      image.read(blob);
#else
      ostringstream pname;
      pname << input_file << "[" << l << "]";
#pragma omp critical
      {
	image.read(pname.str());
      }
#endif
      image.modifyImage();
      bool adaptive = convert_to_gray(image, invert, adaptive_option, verbose);

    
      
      vector<vector<string> > array_of_structures(num_resolutions);
      vector<vector<double> > array_of_avg_bonds(num_resolutions), array_of_ind_conf(num_resolutions);
      vector<vector<Image> > array_of_images(num_resolutions);
      vector<vector<box_t> > array_of_boxes(num_resolutions);

     

      if (input_resolution > 300)
        {
          int percent = (100 * 300) / input_resolution;
          ostringstream scale;
          scale << percent << "%";
          image.scale(scale.str());
          page_scale /= (double) percent / 100;
        }

      if (verbose)
        cout << "Input resolutions are " << select_resolution << endl;

      ColorGray bgColor = getBgColor(image);
      if (rotate != 0)
        {
          image.backgroundColor(bgColor);
          image.rotate(rotate);
        }

      double rotation = 0;
      int unpaper_dx = 0;
      int unpaper_dy = 0;
      for (int i = 0; i < do_unpaper; i++)
        {
          double radians=0;
          int dx=0, dy=0;
          unpaper(image,radians,dx,dy);
          rotation +=radians;
          unpaper_dx +=dx;
          unpaper_dy +=dy;
        }

      // 0.1 is used for THRESHOLD_BOND here to allow for farther processing.
      list<list<list<point_t> > > clusters = find_segments(image, 0.1, bgColor, adaptive, is_reaction, arrows[l], pluses[l], verbose);

      if (verbose)
        cout << "Number of clusters: " << clusters.size() << '.' << endl;

      vector<box_t> boxes;
      int n_boxes = prune_clusters(clusters, boxes);
      std::sort(boxes.begin(), boxes.end(), comp_boxes);

      if (verbose)
        cout << "Number of boxes: " << boxes.size() << '.' << endl;

      // This will hide the output "Warning: non-positive median line gap" from GOCR. Remove after this is fixed:
      fclose(stderr);
      OpenBabel::obErrorLog.StopLogging();

      for (int res_iter = 0; res_iter < num_resolutions; res_iter++)
        {
          int total_boxes = 0;
          double total_confidence = 0;

          int resolution = select_resolution[res_iter];
          int working_resolution = resolution;
          if (resolution > 300)
            working_resolution = 300;

          double THRESHOLD_BOND = set_threshold(threshold,resolution);

          int max_font_height = MAX_FONT_HEIGHT * working_resolution / 150;
          int max_font_width = MAX_FONT_WIDTH * working_resolution / 150;
          bool thick = true;
          if (resolution < 150)
            thick = false;
          else if (resolution == 150 && !jaggy)
            thick = false;


          //Image dbg = image;
          //dbg.modifyImage();
          //dbg.backgroundColor("white");
          //dbg.erase();
          //dbg.type(TrueColorType);
          for (int k = 0; k < n_boxes; k++)
            if ((boxes[k].x2 - boxes[k].x1) > max_font_width && (boxes[k].y2 - boxes[k].y1) > max_font_height
                && !boxes[k].c.empty() && ((boxes[k].x2 - boxes[k].x1) > 2 * max_font_width || (boxes[k].y2
                                           - boxes[k].y1) > 2 * max_font_height))
              {
                int n_atom = 0, n_bond = 0, n_letters = 0, n_label = 0;
                vector<atom_t> atom;
                vector<bond_t> bond;
                vector<letters_t> letters;
                vector<label_t> label;
                double box_scale = 1;
                Image orig_box(Geometry(boxes[k].x2 - boxes[k].x1 + 2 * FRAME, boxes[k].y2 - boxes[k].y1 + 2
                                        * FRAME), bgColor);

                for (unsigned int p = 0; p < boxes[k].c.size(); p++)
                  {
                    int x = boxes[k].c[p].x;
                    int y = boxes[k].c[p].y;
                    ColorGray color = image.pixelColor(x, y);
                    //dbg.pixelColor(x, y, color);
                    orig_box.pixelColor(x - boxes[k].x1 + FRAME, y - boxes[k].y1 + FRAME, color);
                  }

                int width = orig_box.columns();
                int height = orig_box.rows();
                Image thick_box;
                create_thick_box(orig_box,thick_box,width,height,resolution,working_resolution,box_scale,bgColor,THRESHOLD_BOND,res_iter,thick,jaggy);

                if (verbose)
                  cout << "Analysing box " << boxes[k].x1 << "x" << boxes[k].y1 << "-" << boxes[k].x2 << "x" << boxes[k].y2 << " using working resolution " << working_resolution << '.' << endl;

                Image box;
                if (thick)
                  box = thin_image(thick_box, THRESHOLD_BOND, bgColor);
                else
                  box = thick_box;
                potrace_state_t * const  st = raster_to_vector(box,bgColor,THRESHOLD_BOND,width,height,working_resolution);
                potrace_path_t const * const p = st->plist;
                n_atom = find_atoms(p, atom, bond, &n_bond,width,height);

                int real_font_width, real_font_height;
                n_letters = find_chars(p, orig_box, letters, atom, bond, n_atom, n_bond, height, width, bgColor,
                                       THRESHOLD_BOND, max_font_width, max_font_height, real_font_width, real_font_height,verbose);
                if (verbose)
                  cout << "Number of atoms: " << n_atom << ", bonds: " << n_bond << ", " << n_letters << " letters: " << n_letters << " " << letters << " after find_atoms()" << endl;

                double avg_bond_length = percentile75(bond, n_bond, atom);

                double max_area = avg_bond_length * 5;
                if (thick)
                  max_area = avg_bond_length;
                n_letters = find_plus_minus(p, orig_box, bgColor, THRESHOLD_BOND, letters, atom, bond, n_atom, n_bond, height, width,
                                            real_font_height, real_font_width, n_letters, avg_bond_length);
                n_atom = find_small_bonds(p, atom, bond, n_atom, &n_bond, max_area, avg_bond_length / 2, 5);

		//remove_small_bonds_in_chars(atom,bond,letters);


                find_old_aromatic_bonds(p, bond, n_bond, atom, n_atom, avg_bond_length);

                if (verbose)
                  cout << "Number of atoms: " << n_atom << ", bonds: " << n_bond << ", " << n_letters << "letters: " << letters << " after find_old_aromatic_bonds()" << endl;

                double dist = 3.;
                if (working_resolution < 150)
                  dist = 2;

                double thickness = skeletize(atom, bond, n_bond, box, THRESHOLD_BOND, bgColor, dist, avg_bond_length);
                remove_disconnected_atoms(atom, bond, n_atom, n_bond);
                collapse_atoms(atom, bond, n_atom, n_bond, 3);
                remove_zero_bonds(bond, n_bond, atom);

		n_bond = find_wavy_bonds(bond,n_bond,atom,avg_bond_length);
		//				if (ttt++ == 0)  debug_image(orig_box, atom, n_atom, bond, n_bond, "tmp.png");                
                n_letters = find_fused_chars(bond, n_bond, atom, letters, n_letters, real_font_height,
                                             real_font_width, 0, orig_box, bgColor, THRESHOLD_BOND, 3, verbose);

                n_letters = find_fused_chars(bond, n_bond, atom, letters, n_letters, real_font_height,
                                             real_font_width, '*', orig_box, bgColor, THRESHOLD_BOND, 5, verbose);

                flatten_bonds(bond, n_bond, atom, 3);
                remove_zero_bonds(bond, n_bond, atom);
                avg_bond_length = percentile75(bond, n_bond, atom);

                if (verbose)
                  cout << "Average bond length: " << avg_bond_length << endl;

                double max_dist_double_bond = dist_double_bonds(atom, bond, n_bond, avg_bond_length);
                n_bond = double_triple_bonds(atom, bond, n_bond, avg_bond_length, n_atom, max_dist_double_bond);
                n_atom = find_dashed_bonds(p, atom, bond, n_atom, &n_bond, max(MAX_DASH, int(avg_bond_length / 3)),
                                           avg_bond_length, orig_box, bgColor, THRESHOLD_BOND, thick, avg_bond_length, letters);

		

                n_letters = remove_small_bonds(bond, n_bond, atom, letters, n_letters, real_font_height,
                                               MIN_FONT_HEIGHT, avg_bond_length);
		
		n_letters = find_numbers(p, orig_box, letters, atom, bond, n_atom, n_bond, height, width, bgColor,
					 THRESHOLD_BOND, n_letters);

                dist = 4.;
                if (working_resolution < 300)
                  dist = 3;
                if (working_resolution < 150)
                  dist = 2;

                n_bond = fix_one_sided_bonds(bond, n_bond, atom, dist, avg_bond_length);

                n_letters = clean_unrecognized_characters(bond, n_bond, atom, real_font_height, real_font_width, 4,
		          letters, n_letters);

                thickness = find_wedge_bonds(thick_box, atom, n_atom, bond, n_bond, bgColor, THRESHOLD_BOND,
                                             max_dist_double_bond, avg_bond_length, 3, 1);

                n_label = assemble_labels(letters, n_letters, label);
		

                if (verbose)
                  cout << n_label << " labels: " << label << " after assemble_labels()" << endl;

                remove_disconnected_atoms(atom, bond, n_atom, n_bond);
	

                collapse_atoms(atom, bond, n_atom, n_bond, thickness);

                remove_zero_bonds(bond, n_bond, atom);

                flatten_bonds(bond, n_bond, atom, 2 * thickness);

                remove_zero_bonds(bond, n_bond, atom);

                avg_bond_length = percentile75(bond, n_bond, atom);

                collapse_double_bonds(bond, n_bond, atom, max_dist_double_bond);

	

                extend_terminal_bond_to_label(atom, letters, n_letters, bond, n_bond, label, n_label, avg_bond_length / 2,
					      thickness, max_dist_double_bond);

		

                remove_disconnected_atoms(atom, bond, n_atom, n_bond);
                collapse_atoms(atom, bond, n_atom, n_bond, thickness);
                collapse_doubleup_bonds(bond, n_bond);

	
                remove_zero_bonds(bond, n_bond, atom);
                flatten_bonds(bond, n_bond, atom, thickness);
                remove_zero_bonds(bond, n_bond, atom);
                remove_disconnected_atoms(atom, bond, n_atom, n_bond);

	
                
                extend_terminal_bond_to_bonds(atom, bond, n_bond, avg_bond_length, 2 * thickness, max_dist_double_bond);




                collapse_atoms(atom, bond, n_atom, n_bond, 3);


                remove_zero_bonds(bond, n_bond, atom);
                flatten_bonds(bond, n_bond, atom, 5);
                remove_zero_bonds(bond, n_bond, atom);
	
                n_letters = clean_unrecognized_characters(bond, n_bond, atom, real_font_height, real_font_width, 0,
                            letters, n_letters);
		int recognized_chars = count_recognized_chars(atom,bond);
	
                assign_charge(atom, bond, n_atom, n_bond, spelling, superatom, debug);
                find_up_down_bonds(bond, n_bond, atom, thickness);
                int real_atoms = count_atoms(atom, n_atom);
                int bond_max_type = 0;
                int real_bonds = count_bonds(bond, n_bond,bond_max_type);

                if (verbose)
                  cout << "Final number of atoms: " << real_atoms << ", bonds: " << real_bonds << ", chars: " << n_letters << '.' << endl;

                split_fragments_and_assemble_structure_record(atom,n_atom,bond,n_bond,boxes,
							      l,k,resolution,res_iter,output_image_file_prefix,image,orig_box,real_font_width,real_font_height,
							      thickness,avg_bond_length,superatom,real_atoms,real_bonds,bond_max_type,
							      box_scale,page_scale,rotation,unpaper_dx,unpaper_dy,output_format,embedded_format,is_reaction,show_confidence,
							      show_resolution_guess,show_page,show_coordinates, show_avg_bond_length,array_of_structures,
							      array_of_avg_bonds,array_of_ind_conf,array_of_images,array_of_boxes,total_boxes,total_confidence,
							      recognized_chars,show_learning,res_iter,verbose);

                if (st != NULL)
                  potrace_state_free(st);
              }
	  array_of_confidence[res_iter] += total_confidence;
	  boxes_per_res[res_iter] += total_boxes;
          //dbg.write("debug.png");
        }

     
      
      #pragma omp critical
      {
         if (show_learning)
	  for (int j = 0; j < num_resolutions; j++)
	    for (unsigned int i = 0; i < array_of_structures[j].size(); i++)
	    {
	      pages_of_structures[l].push_back(array_of_structures[j][i]);
	      if (!output_image_file_prefix.empty())
		pages_of_images[l].push_back(array_of_images[j][i]);
	      pages_of_avg_bonds[l].push_back(array_of_avg_bonds[j][i]);
	      pages_of_ind_conf[l].push_back(array_of_ind_conf[j][i]);
	      pages_of_boxes[l].push_back(array_of_boxes[j][i]);
	      total_structure_count++;
	    }
	 else
	   for (int j = 0; j < num_resolutions; j++)
	     {
                array_of_structures_page[l][j] = array_of_structures[j];
		if (!output_image_file_prefix.empty())
		  array_of_images_page[l][j] = array_of_images[j];
		array_of_avg_bonds_page[l][j] = array_of_avg_bonds[j];
		array_of_ind_conf_page[l][j] = array_of_ind_conf[j];
		array_of_boxes_page[l][j] = array_of_boxes[j];
	      }
	   
       }
     }

    double max_conf = -FLT_MAX;
    int max_res = 0;
    for (int i = 0; i < num_resolutions; i++)
        {
          if (boxes_per_res[i] > 0 && array_of_confidence[i]/boxes_per_res[i] > max_conf)
            {
              max_conf = array_of_confidence[i]/boxes_per_res[i];
              max_res = i;
            }
        }
      for (int i = 0; i < num_resolutions; i++)
	if (boxes_per_res[i] > 0 && array_of_confidence[i]/boxes_per_res[i] == max_conf && select_resolution[i] == 300) // second 300 dpi is without thinning
	  {
	    max_res = i;
	    break;
	  }

	if (!show_learning)
         for (int l = 0; l < page; l++) 
	    {
	      pages_of_structures[l] = array_of_structures_page[l][max_res];
	      if (!output_image_file_prefix.empty())
		pages_of_images[l] = array_of_images_page[l][max_res];
	      pages_of_avg_bonds[l] = array_of_avg_bonds_page[l][max_res];
	      pages_of_ind_conf[l] = array_of_ind_conf_page[l][max_res];
	      pages_of_boxes[l] = array_of_boxes_page[l][max_res];
	      total_structure_count += array_of_structures_page[l][max_res].size();
	    }

  double best_bond = 0;

  //if (total_structure_count >= STRUCTURE_COUNT)
  //  find_limits_on_avg_bond(best_bond, pages_of_avg_bonds, pages_of_ind_conf);

  // If multiple pages are processed at several  resolutions different pages
  // may be processed at different resolutions leading to a seemingly different average bond length
  // Currently multi-page documents (PDF and PS) are all processed at the same resolution
  // and single-page images have all structures on the page at the same resolution

  //cout << min_bond << " " << max_bond << endl;

#ifdef OSRA_LIB
  ostream &out_stream = structure_output_stream;
#else
  ostream &out_stream = outfile.is_open() ? outfile : cout;
#endif


  // For Andriod version we will find the structure with maximum confidence value, as the common usecase for Andriod is to analyse the
  // image (taken by embedded photo camera) that usually contains just one molecule:
  double max_confidence = -FLT_MAX;
  int l_index = 0;
  int i_index = 0;


  int image_count = 0;

  for (int l = 0; l < page; l++)
    {
      for (unsigned int i = 0; i < pages_of_structures[l].size(); i++)
	if (best_bond == 0 || (pages_of_avg_bonds[l][i] > best_bond/2 && pages_of_avg_bonds[l][i] < 2*best_bond))
	  {
	    if (pages_of_ind_conf[l][i] > max_confidence)
	      {
		max_confidence = pages_of_ind_conf[l][i];
		l_index = l;
		i_index = i;
	      }
	    
	    if (output_format != "mol" && !is_reaction)
	      {
		out_stream << pages_of_structures[l][i];

		// Dump this structure into a separate file:
		if (!output_image_file_prefix.empty())
		  {
		    ostringstream fname;
		    fname << output_image_file_prefix << image_count << ".png";
		    image_count++;
		    if (fname.str() != "")
		      {
			Image tmp = pages_of_images[l][i];
			if (resize != "")
			  {
			    tmp.scale(resize);
			  }
			tmp.write(fname.str());
		      }
		  }
	      }
	  }
      if (is_reaction && !arrows[l].empty())
	 {
	   vector<string> reactions;
	   vector<box_t> rbox;
	   arrange_reactions(arrows[l], pages_of_boxes[l], pluses[l], reactions, rbox, pages_of_structures[l],output_format);
	   for (int k=0; k<reactions.size(); k++)
	     {
	       out_stream << reactions[k]<<endl;

	       if (!output_image_file_prefix.empty())
		 {
		   ostringstream fname;
		   fname << output_image_file_prefix << image_count << ".png";
		   image_count++;
		   if (fname.str() != "")
		     {
		       Image tmp = pages_of_images[l][k];
		       Geometry geometry = Geometry(rbox[k].x2 - rbox[k].x1, rbox[k].y2 - rbox[k].y1, rbox[k].x1, rbox[k].y1);
		       tmp.crop(geometry);
		       if (resize != "")
			 {
			   tmp.scale(resize);
			 }
		       tmp.write(fname.str());
		     }
		 }
	     }
	 }
    }
  // Output the structure with maximum confidence value:
  if (output_format == "mol")
    {
      out_stream << pages_of_structures[l_index][i_index];
      if (!output_image_file_prefix.empty())
	{
	  ostringstream fname;
	  fname << output_image_file_prefix  << ".png";
	  if (fname.str() != "")
	    {
	      Image tmp = pages_of_images[l_index][i_index];
	      if (resize != "")
		{
		  tmp.scale(resize);
		}
	      tmp.write(fname.str());
	    }
	}
    }


  out_stream.flush();

#ifndef OSRA_LIB
  if (!output_file.empty())
    outfile.close();
#endif

  return 0;
}
