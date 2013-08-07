package cadd.osra.main;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URLEncoder;

import android.app.Activity;
import android.content.Intent;
import android.net.Uri;
import android.os.Bundle;
import android.widget.TextView;

public class runosra extends Activity {
    /** Called when the activity is first created. */
    @Override
    public void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.main);
 
        try {
        	writeToStream(getAssets().open("spelling.txt"), openFileOutput("spelling.txt", MODE_WORLD_READABLE ));
        } catch (Exception e) {
			e.printStackTrace();
		}
        try {
        	writeToStream(getAssets().open("superatom.txt"), openFileOutput("superatom.txt", MODE_WORLD_READABLE));
        } catch (Exception e) {
			e.printStackTrace();
		}
        byte [] rawimage=null;
        try {
        	rawimage=writeToArray(getAssets().open("c.jpg"));
        } catch (Exception e) {
			e.printStackTrace();
		}
        
        
            
    	File spelling=getFileStreamPath("spelling.txt");
    	File superatom=getFileStreamPath("superatom.txt");
    	//File image=getFileStreamPath("chemnav.png");
        //TextView  tv = new TextView(this);
        String [] jargv= {"osra","-f","inchi","-l",spelling.getAbsolutePath(),"-a",superatom.getAbsolutePath()};
        String inchi=nativeosra(jargv,rawimage);
        String enc_inchi="";
        try {
        	enc_inchi=URLEncoder.encode(inchi, "UTF-8");
        } catch (Exception e) {
			e.printStackTrace();
		}	
        if (enc_inchi.length()!=0) {
        	//String base_url="http://en.wikipedia.org/wiki/Special:Search?fulltext=Search&search=";
        	String base_url="http://129.43.27.140/cgi-bin/lookup/results?type=inchi&context_all=all&query=";
        	String url=base_url+enc_inchi;
        	Uri uri = Uri.parse(url);
        	Intent intent = new Intent(Intent.ACTION_VIEW, uri);
        	startActivity(intent);
        	 
        }
        //tv.setText( enc_inchi );
        //setContentView(tv);
    
    }
   
    public native String  nativeosra(String [] jargv, byte [] rawimage);
    static {
     	System.loadLibrary("openbabel");
        System.loadLibrary("osra");
    }
 
    public static void writeToStream(InputStream in , OutputStream out)	throws IOException 
    {
    	byte[] bytes = new byte[2048];
	
    	for (int c = in.read(bytes); c != -1; c = in.read(bytes)) {
    		out.write(bytes,0, c);
    	}
    	in.close();
    	out.close();    		
    }
    
    public static byte [] writeToArray(InputStream in)	throws IOException 
    {
    	ByteArrayOutputStream out=new ByteArrayOutputStream();
    	byte[] bytes = new byte[2048];
    	
    	for (int c = in.read(bytes); c != -1; c = in.read(bytes)) {
    		out.write(bytes,0, c);
    	}
    	in.close();
    	//byte [] prep = new byte[out.size()];
    	byte [] arr=out.toByteArray();
    	out.close();
    	return arr;
    }
}
