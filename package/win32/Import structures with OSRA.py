dot_version="1.3.5"

import tempfile
import StringIO
import popen2
import ImageGrab
import os
import sys
import ChemScript12

class TempFile(object):
    def __init__(self, *args, **nargs):
        self.des, self.name = tempfile.mkstemp(*args, **nargs)
    def close(self):
        os.close(self.des)
        os.remove(self.name)

class WindowsClip(object):
    def createtempfile(self):
        self.tempfile = TempFile(suffix=".png")
    def getimage(self):
        self.image = ImageGrab.grabclipboard()
        return self.image is not None
    def saveimage(self):
        self.image.save(self.tempfile.name)
    def deletetempfile(self):
        self.tempfile.close()

def osra(filename):
    osra = os.environ.get("OSRA", None)
    osra = os.path.join(osra,'osra.bat')
    if not os.path.isfile(osra):
        osra = os.path.join(os.path.realpath(os.path.dirname(sys.argv[0])),'osra.bat')
    if not os.path.isfile(osra):
        programfiles = os.environ.get("PROGRAMFILES", None)
        osra = os.path.join(programfiles, "osra", dot_version,"osra.bat")
    stdout, stdin, stderr = popen2.popen3('"%s" -f sdf %s' % (osra, filename))  
    sdf = stdout.read()
    return sdf
	
clip = WindowsClip()
success = clip.getimage()
if success:  
    clip.createtempfile()
    clip.saveimage()
    sdf = osra(clip.tempfile.name)
    clip.deletetempfile()
    if  sdf.endswith("$$$$\n"):
        m1 = ChemScript12.StructureData.LoadData(sdf)
        m1.WriteFile(sys.argv[2])



