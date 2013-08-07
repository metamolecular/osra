"""Authors: Noel O'Boyle and Igor V. Filippov
Copied of......hmmm... Inspired by the "fetch from webbook" plugin :-)
Converts an image from clipboard to a molecule using OSRA    
"""

import os
import popen2
import oasa_bridge
import dialogs
import tempfile
import Pmw
import StringIO
import os, sys


def err_mess_box(mess, title="Error"): #Pops up error OK-box
	message = ""
	for m in mess:
		message=message+m+"\n"
	dialog = Pmw.Dialog(App.paper, buttons=('OK',),
	defaultbutton='OK', title=title)

	w = Pmw.LabeledWidget(dialog.interior(), labelpos='n', label_text=message)
	w.pack(expand=1, fill='both', padx=4, pady=4)
	dialog.activate()
	
def run_osra(osra):
	sdf = " "    
	filedes, filename = tempfile.mkstemp(suffix='.png')

	if os.name=="posix":
		import pygtk
		pygtk.require('2.0')
		import gtk, gobject
		clipboard = gtk.clipboard_get()
		image=clipboard.wait_for_image()
		if not image:
			return sdf
		try:
			image.save(filename,"png")
		except:
			return sdf
	else:
		import ImageGrab
		image = ImageGrab.grabclipboard()
		if not image:
			return sdf
		try:
			image.save(filename)
		except:
			return sdf

	try:
		stdout, stdin, stderr = popen2.popen3('"%s" -f sdf %s' % (osra, filename))
	except:
		os.remove(filename)
		return sdf

	sdf = stdout.read()
	#os.remove(filename)
	return sdf


def present_mol(sdf):
	if not sdf.rstrip().endswith("$$$$"):
		return 0
	try:
		mol = StringIO.StringIO(sdf)
		molec = oasa_bridge.read_molfile(mol, App.paper)
		mol.close()
	except:
		return 0
	if len(molec.atoms)<2:
		return 0
	averagey = sum([atom.y for atom in molec.atoms]) / float(len(molec.atoms))
	for atom in molec.atoms:
		atom.y =  2*averagey - atom.y
	N = 0
	for minimol in molec.get_disconnected_subgraphs():
		N += 1
		App.paper.stack.append(minimol)
		minimol.draw()
		App.paper.add_bindings()
		App.paper.start_new_undo_record()
	return N



osra = os.environ.get("OSRA", None)
if osra and  os.path.isfile(osra):
	if os.name=="posix":
		r, w = os.pipe() # these are file descriptors, not file objects
		pid = os.fork()
		if pid:
			# we are the parent
			os.close(w) # use os.close() to close a file descriptor
			r = os.fdopen(r) # turn r into a file object
			dialog = dialogs.progress_dialog(App, title="Progress")
			dialog.update(0.1, top_text = "Calling OSRA...", bottom_text = "Image processing in progress")
			sdf = r.read()
			dialog.update(0.9, top_text = "Adding molecules to workspace...", bottom_text = "Almost there!")
			N = present_mol(sdf)
			dialog.close()
			if N<1:
				err_mess_box(["Image could not be converted to a molecule."])
#			else:
#				err_mess_box(["%d molecule%s added" % (N, ["s", ""][N==1])], "Info")
			os.waitpid(pid, 0) # make sure the child process gets cleaned up
		else:
			# we are the child
			os.close(r)
			w = os.fdopen(w, 'w')
			sdf = run_osra(osra)
			w.write(sdf)
			w.close()
			sys.exit(0)
	else:
		dialog = dialogs.progress_dialog(App, title="Progress")
		dialog.update(0.1, top_text = "Calling OSRA...", bottom_text = "Image processing in progress")
		sdf = run_osra(osra)
		dialog.update(0.9, top_text = "Adding molecules to workspace...", bottom_text = "Almost there!")
		N = present_mol(sdf)
		dialog.close()
		if N<1:
			err_mess_box(["Image could not be converted to a molecule."])
#		else:
#			err_mess_box(["%d molecule%s added" % (N, ["s", ""][N==1])], "Info")
		
else:
    err_mess_box(["You need to set the environment variable " \
		  "OSRA to point to the OSRA executable.\n" \
                  "When setting the variable, do not include quotation " \
                  "marks around the path."])


