package net.sf.osra;

import java.io.Writer;

/**
 * JNI bridge for OSRA library.
 * 
 */

public class OsraLib {
	
	/**
	 * Process the given image with OSRA library. For more information see the corresponding <a href=
	 * "https://sourceforge.net/apps/mediawiki/osra/index.php?title=Usage">CLI options</a>.
	 * 
	 * @param imageData
	 *            the image binary data
	 * @param outputStructureWriter
	 *            the writer to output the found structures in given format
	 * @param rotate
	 *            rotate image, degrees
	 * @param invert
	 *            force color inversion (for white-on-black images)
	 * @param input_resolution
	 *             force processing at a specific resolution, dpi
	 * @param threshold
	 *              black-white binarization threshold 0.0-1.0
	 * @param do_unpaper
	 *              perform unpaper image pre-processing, rounds
	 * @param jaggy
	 *              perform image downsampling
	 * @param adaptive_option
	 *              perform adaptive thresholding (more CPU-intensive)
	 * @param format
	 *            one of the formats, accepted by OpenBabel ("sdf", "smi", "can").
	 * @param embeddedFormat
	 *            format to be embedded into SDF ("inchi", "smi", "can").
	 * @param outputConfidence
	 *            include confidence
	 * @param show_resolution_guess
	 *            include image resolution estimate
	 * @param show_page
	 *            include page number
	 * @param outputCoordinates
	 *            include box coordinates
	 * @param outputAvgBondLength
	 *            include average bond length
	 * @return 0, if the call succeeded or negative value in case of error
	 */
    public static native int processImage(byte[] imageData, Writer outputStructureWriter,   
					  int rotate, boolean invert, int input_resolution, double threshold, int do_unpaper, boolean jaggy, boolean adaptive_option,
					  String format,
					  String embeddedFormat, 
					  boolean outputConfidence,
					  boolean show_resolution_guess,
					  boolean show_page,
					  boolean outputCoordinates, 
					  boolean outputAvgBondLength);

	private static final String NAME = "osra_java";
		
	static {
		try {
			System.out.println("Osra Name2: "+NAME);
			System.loadLibrary(NAME);
		} catch (Exception e) {
			throw new RuntimeException(e);
		} 
	}

}

