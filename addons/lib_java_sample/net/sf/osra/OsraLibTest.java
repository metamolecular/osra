package net.sf.osra;

/*import static org.junit.Assert.assertNotNull;*/

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.StringWriter;
import java.io.InputStream;
import java.io.BufferedInputStream;
import java.io.ByteArrayOutputStream;

/*import net.sf.osra.OsraLib;*/

/*import org.apache.commons.io.IOUtils;*/
/*import org.junit.Test;*/

/**
 * Sample usage of the library.
 * 
 * @author <a href="mailto:dmitry.katsubo@gmail.com">Dmitry Katsubo</a>
 */
public class OsraLibTest {



	public void testProcessImage() throws IOException {
		final StringWriter writer = new StringWriter();
		InputStream is = new BufferedInputStream(new FileInputStream("test.png"));
		byte[] imageData = writeToArray(is);
		int result = OsraLib.processImage(imageData,
				writer, 0, false, 0, 0, 0, 
				false, false,"sdf", "inchi", true, true, true, true, true);

		System.out.println("OSRA completed with result:" + result + " structure:\n" + writer.toString() + "\n");
	}

	public static void main(String[] args) throws IOException {
		new OsraLibTest().testProcessImage();
	}


 public static byte [] writeToArray(InputStream in)	throws IOException  {
    	ByteArrayOutputStream out=new ByteArrayOutputStream();
    	byte[] bytes = new byte[2048];
    	
    	for (int c = in.read(bytes); c != -1; c = in.read(bytes)) {
    		out.write(bytes,0, c);
    	}
    	in.close();
    	byte [] arr=out.toByteArray();
    	out.close();
    	return arr;
    }
}
