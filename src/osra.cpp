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

#include <string.h> // strncpy()
#include <libgen.h> // dirname()

#include <tclap/CmdLine.h>

#include "osra_lib.h"
#include "config.h" // PACKAGE_VERSION

int main(int argc,
         char **argv
        )
{
  TCLAP::CmdLine cmd("OSRA: Optical Structure Recognition Application, created by Igor Filippov, 2013", ' ',
                     PACKAGE_VERSION);

  //
  // Image pre-processing options
  //
  TCLAP::ValueArg<double> rotate_option("R", "rotate", "Rotate image clockwise by specified number of degrees", false, 0,
                                        "0..360");
  cmd.add(rotate_option);

  TCLAP::SwitchArg invert_option("n", "negate", "Invert color (white on black)", false);
  cmd.add(invert_option);

  TCLAP::ValueArg<int> resolution_option("r", "resolution", "Resolution in dots per inch", false, 0, "default: auto");
  cmd.add(resolution_option);

  TCLAP::ValueArg<double> threshold_option("t", "threshold", "Gray level threshold", false, 0, "0.2..0.8");
  cmd.add(threshold_option);

  TCLAP::ValueArg<int> do_unpaper_option("u", "unpaper", "Pre-process image with unpaper algorithm, rounds", false, 0,
                                         "default: 0 rounds");
  cmd.add(do_unpaper_option);

  TCLAP::SwitchArg jaggy_option("j", "jaggy", "Additional thinning/scaling down of low quality documents", false);
  cmd.add(jaggy_option);

  TCLAP::SwitchArg adaptive_option("i", "adaptive", "Adaptive thresholding pre-processing, useful for low light/low contrast images", false);
  cmd.add(adaptive_option);

  //
  // Output format options
  //
  TCLAP::ValueArg<string> output_format_option("f", "format", "Output format", false, "can", "can/smi/sdf");
  cmd.add(output_format_option);

  TCLAP::ValueArg<string> embedded_format_option("", "embedded-format", "Embedded format", false, "", "inchi/smi/can");
  cmd.add(embedded_format_option);

  TCLAP::SwitchArg show_confidence_option("p", "print", "Print out confidence estimate", false);
  cmd.add(show_confidence_option);

  TCLAP::SwitchArg show_resolution_guess_option("g", "guess", "Print out resolution guess", false);
  cmd.add(show_resolution_guess_option);

  TCLAP::SwitchArg show_page_option("e", "page", "Show page number for PDF/PS/TIFF documents (only for SDF/SMI/CAN output format)", false);
  cmd.add(show_page_option);

  TCLAP::SwitchArg show_coordinates_option("c", "coordinates", "Show surrounding box coordinates (only for SDF/SMI/CAN output format)", false);
  cmd.add(show_coordinates_option);

  TCLAP::SwitchArg show_avg_bond_length_option("b", "bond", "Show average bond length in pixels (only for SDF/SMI/CAN output format)", false);
  cmd.add(show_avg_bond_length_option);

  //
  // Dictionaries options
  //
  TCLAP::ValueArg<string> spelling_file_option("l", "spelling", "Spelling correction dictionary", false, "", "configfile");
  cmd.add(spelling_file_option);

  TCLAP::ValueArg<string> superatom_file_option("a", "superatom", "Superatom label map to SMILES", false, "", "configfile");
  cmd.add(superatom_file_option);

  //
  // Debugging options
  //
  TCLAP::SwitchArg debug_option("d", "debug", "Print out debug information on spelling corrections", false);
  cmd.add(debug_option);

  TCLAP::SwitchArg verbose_option("v", "verbose", "Be verbose and print the program flow", false);
  cmd.add(verbose_option);

  TCLAP::ValueArg<string> output_image_file_prefix_option("o", "output", "Write recognized structures to image files with given prefix", false, "", "filename prefix");
  cmd.add(output_image_file_prefix_option);

  TCLAP::ValueArg<string> resize_option("s", "size", "Resize image on output", false, "", "dimensions, 300x400");
  cmd.add(resize_option);

  //
  // Input-output options
  //
  TCLAP::UnlabeledValueArg<string> input_file_option("in", "input file", true, "", "filename");
  cmd.add(input_file_option);

  TCLAP::ValueArg<string> output_file_option("w", "write", "Write recognized structures to text file", false, "", "filename");
  cmd.add(output_file_option);

  TCLAP::SwitchArg show_learning_option("", "learn", "Print out all structure guesses with confidence parameters", false);
  cmd.add(show_learning_option);

  cmd.parse(argc, argv);

  // Calculating the current dir:
  char progname[1024];
  strncpy(progname, cmd.getProgramName().c_str(), sizeof(progname) - 1);
  progname[sizeof(progname) - 1] = '\0';
  string osra_dir = dirname(progname);

  int result = osra_process_image(
                 input_file_option.getValue(),
                 output_file_option.getValue(),
                 rotate_option.getValue(),
                 invert_option.getValue(),
                 resolution_option.getValue(),
                 threshold_option.getValue(),
                 do_unpaper_option.getValue(),
                 jaggy_option.getValue(),
                 adaptive_option.getValue(),
                 output_format_option.getValue(),
                 embedded_format_option.getValue(),
                 show_confidence_option.getValue(),
                 show_resolution_guess_option.getValue(),
                 show_page_option.getValue(),
                 show_coordinates_option.getValue(),
                 show_avg_bond_length_option.getValue(),
		 show_learning_option.getValue(),
                 osra_dir,
                 spelling_file_option.getValue(),
                 superatom_file_option.getValue(),
                 debug_option.getValue(),
                 verbose_option.getValue(),
                 output_image_file_prefix_option.getValue(),
                 resize_option.getValue()
               );

  return result;
}
