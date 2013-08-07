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

/* Fix for jlong definition in jni.h on some versions of gcc on Windows */
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
typedef long long __int64;
#endif

#include <jni.h>

#include <stdlib.h> // calloc(), free()

#include <string> // std::string
#include <ostream> // std:ostream
#include <sstream> // std:ostringstream

#include "config.h" // PACKAGE_VERSION

using namespace std;


#ifdef OSRA_JAVA
#include "osra_lib.h"

extern "C" {
  /*
   * Class:     net_sf_osra_OsraLib
   * Method:    processImage
   * Signature: ([BLjava/io/Writer;Ljava/lang/String;Ljava/lang/String;ZZZ)I
   */
  JNIEXPORT jint JNICALL Java_net_sf_osra_OsraLib_processImage(JNIEnv *, jclass, jbyteArray, jobject, jint, jboolean,jint,jdouble,jint, jboolean, jboolean,jstring, jstring,  jboolean, jboolean,jboolean, jboolean, jboolean);

  /*
   * Class:     net_sf_osra_OsraLib
   * Method:    getVersion
   * Signature: ()Ljava/lang/String;
   */
  JNIEXPORT jstring JNICALL Java_net_sf_osra_OsraLib_getVersion(JNIEnv *, jclass);
}

JNIEXPORT jint JNICALL Java_net_sf_osra_OsraLib_processImage(JNIEnv *j_env, jclass j_class,
    jbyteArray j_image_data,
    jobject j_writer,
    jint j_rotate,
    jboolean j_invert,
    jint j_input_resolution,
    jdouble j_threshold,
    jint j_do_unpaper,
    jboolean j_jaggy,
    jboolean j_adaptive_option,
    jstring j_output_format,
    jstring j_embedded_format,
    jboolean j_output_confidence,
    jboolean j_show_resolution_guess,
    jboolean j_show_page,
    jboolean j_output_coordinates,
    jboolean j_output_avg_bond_length)
{
  const char *output_format = j_env->GetStringUTFChars(j_output_format, NULL);
  const char *embedded_format = j_env->GetStringUTFChars(j_embedded_format, NULL);
  const char *image_data = (char *) j_env->GetByteArrayElements(j_image_data, NULL);

  int result = -1;

  if (image_data != NULL)
    {
      // Perhaps there is a more optimal way to bridge from std:ostream to java.io.Writer.
      // See http://stackoverflow.com/questions/524524/creating-an-ostream/524590#524590
      ostringstream structure_output_stream;

      result = osra_process_image(
                 image_data,
                 j_env->GetArrayLength(j_image_data),
                 structure_output_stream,
                 j_rotate,
                 j_invert,
                 j_input_resolution,
                 j_threshold,
                 j_do_unpaper,
                 j_jaggy,
                 j_adaptive_option,
                 output_format,
                 embedded_format,
                 j_output_confidence,
                 j_show_resolution_guess,
                 j_show_page,
                 j_output_coordinates,
                 j_output_avg_bond_length,
                 "."
               );

      j_env->ReleaseByteArrayElements(j_image_data, (jbyte *) image_data, JNI_ABORT);

      // Locate java.io.Writer#write(String) method:
      jclass j_writer_class = j_env->FindClass("java/io/Writer");
      jmethodID write_method_id = j_env->GetMethodID(j_writer_class, "write", "(Ljava/lang/String;)V");

      jstring j_string = j_env->NewStringUTF(structure_output_stream.str().c_str());

      j_env->CallVoidMethod(j_writer, write_method_id, j_string);

      j_env->DeleteLocalRef(j_writer_class);
      j_env->DeleteLocalRef(j_string);
    }

  j_env->ReleaseStringUTFChars(j_output_format, output_format);
  j_env->ReleaseStringUTFChars(j_embedded_format, embedded_format);

  return result;
}

JNIEXPORT jstring JNICALL Java_net_sf_osra_OsraLib_getVersion(JNIEnv *j_env, jclass j_class)
{
  return j_env->NewStringUTF(PACKAGE_VERSION);
}
#endif
