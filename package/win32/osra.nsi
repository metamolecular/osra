!define DOT_VERSION  "2.0.0"
!define DASH_VERSION "2-0-0"

!define PIL_EXE "PIL-1.1.6.win32-py2.5.exe"
!define PIL_URL "http://effbot.org/media/downloads/${PIL_EXE}"
!define GS_EXE "gs907w32.exe"
!define GS_URL "http://downloads.ghostscript.com/public/${GS_EXE}"

!include Sections.nsh
; include for some of the windows messages defines
!include "winmessages.nsh"
; HKLM (all users) vs HKCU (current user) defines
!define env_hklm 'HKLM "SYSTEM\CurrentControlSet\Control\Session Manager\Environment"'
!define env_hkcu 'HKCU "Environment"'

!include "TextFunc.nsh" 
!insertmacro LineFind

; The name of the installer
Name "Optical Structure Recognition Application"

; The file to write
OutFile "osra-setup-${DASH_VERSION}.exe"

; The default installation directory
InstallDir $PROGRAMFILES\osra\${DOT_VERSION}

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\osra\${DOT_VERSION}" "Install_Dir"

LicenseData "license.txt"

; Request application privileges for Windows Vista
RequestExecutionLevel admin

;--------------------------------

; Pages

Page license
Page components
Page directory
Page instfiles

UninstPage uninstConfirm
UninstPage instfiles

;--------------------------------

; The stuff to install
Section "osra (required)"

  SectionIn RO
  
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR
  
  ; Put file there
  File "osra-bin.exe"
  File "pthreadGC2.dll"
  File "delegates.mgk"
  File "type.mgk"
  File "type-ghostscript.mgk"
  File "type-solaris.mgk"
  File "type-windows.mgk"
  File "colors.mgk"
  File "log.mgk"
  File "magic.mgk"
  File "modules.mgk"
  File "README.txt"
  File "spelling.txt"
  File "superatom.txt"
  strcpy $3 "GPL Ghostscript"
  call checkSoftVersion
  call getGhostscriptInstPath
  strcmp $1 "" no_gs +2
  no_gs:
  MessageBox MB_OK "Ghostscript interpreter not found" IDOK 0 
  call createOSRAbat
  
  
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\osra\${DOT_VERSION} "Install_Dir" "$INSTDIR"
  WriteRegStr HKLM SOFTWARE\osra "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "DisplayName" "OSRA ${DOT_VERSION}"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra" "NoRepair" 1
  WriteUninstaller "uninstall.exe"
  
  
  ; set variable
  WriteRegExpandStr ${env_hklm} OSRA "$INSTDIR"
  ; make sure windows knows about the change
  SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000

SectionEnd

Section /o "Symyx/Accelrys Draw plugin" symyx_draw
 strcpy $3 "Symyx Technologies, Inc.\Symyx Draw\Client"
 call CheckSoftVersion
 strcmp $2 "" no_symyx
 call getSymyxPath
 strcmp $1 "" no_symyx
 SetOutPath "$1\AddIns"
 File "plugins\symyx_draw\OSRAAction.xml"
 SetOutPath "$1\AddIns\OSRAAction"
 File "plugins\symyx_draw\README.txt"
 File "plugins\symyx_draw\OSRAAction\OSRAAction.dll"
 File "plugins\symyx_draw\OSRAAction\OSRAAction.dll.config"
 Goto done
 no_symyx:
  MessageBox MB_OK "Symyx Draw not found" IDOK done
 done:
SectionEnd

Section /o "ChemBioOffice 12 plugin" chemoffice
; strcpy $3 "CambridgeSoft\ChemScript"
; call CheckSoftVersion
; strcmp $2 "12.0" +1 no_chemoffice
 call getChemScriptPath
 strcmp $1 "" no_chemoffice
 ReadRegStr $2 HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\PIL-py2.5" "DisplayName"
 strcmp $2 "" +1 pil_exists
 call downloadPIL
 pil_exists:
 SetOutPath "$1\Scripts"
 File "plugins\chemoffice\Import Structures with OSRA.py"
 Push "$1\Scripts\Import Structures with OSRA.py" ; file to modify
 Push "dot_version=" ; string that a line must begin with *WS Sensitive*
 Push "dot_version='${DOT_VERSION}'" ; string to replace whole line with
 Call ReplaceLineStr	
 AccessControl::GrantOnFile "$1\Scripts\Import Structures with OSRA.py" "(S-1-5-32-545)" "GenericExecute"
 Goto done
 no_chemoffice:
  MessageBox MB_OK "ChemScript 12.0 not found" IDOK done
 done:
SectionEnd

; Uninstaller

Section "Uninstall"
	# call userInfo plugin to get user info.  The plugin puts the result in the stack
    userInfo::getAccountType
   
    # pop the result from the stack into $0
    pop $0
 
    # compare the result with the string "Admin" to see if the user is admin.
    # If match, jump 3 lines down.
    strCmp $0 "Admin" +3
 
    # if there is not a match, print message and return
    messageBox MB_OK "Please run this with Administrator privileges"
    Quit   
  ReadRegStr $0 HKLM SOFTWARE\osra\${DOT_VERSION} "Install_Dir"
  strcpy $INSTDIR $0
  ; Remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\osra"
  DeleteRegKey HKLM SOFTWARE\osra\${DOT_VERSION}
  ; delete variable
  DeleteRegValue ${env_hklm} OSRA
  ; make sure windows knows about the change
  SendMessage ${HWND_BROADCAST} ${WM_WININICHANGE} 0 "STR:Environment" /TIMEOUT=5000

  ; Remove files and uninstaller
  Delete $INSTDIR\osra-bin.exe
  Delete $INSTDIR\pthreadGC2.dll
  Delete $INSTDIR\delegates.mgk
  Delete $INSTDIR\type.mgk
  Delete $INSTDIR\type-ghostscript.mgk
  Delete $INSTDIR\type-solaris.mgk
  Delete $INSTDIR\type-windows.mgk
  Delete $INSTDIR\colors.mgk
  Delete $INSTDIR\log.mgk
  Delete $INSTDIR\magic.mgk
  Delete $INSTDIR\modules.mgk
  Delete $INSTDIR\README.txt
  Delete $INSTDIR\osra.bat
  Delete $INSTDIR\superatom.txt
  Delete $INSTDIR\spelling.txt
  Delete $INSTDIR\uninstall.exe
  RMDir "$INSTDIR"
  strcpy $3 "Symyx Technologies, Inc.\Symyx Draw\Client"
  call un.CheckSoftVersion
  strcmp $2 "" no_symyx
  call un.getSymyxPath
  strcmp $1 "" no_symyx 
  Delete "$1\AddIns\OSRAAction.xml"
  Delete "$1\AddIns\OSRAAction\README.txt"
  Delete "$1\AddIns\OSRAAction\OSRAAction.dll"
  Delete "$1\AddIns\OSRAAction\OSRAAction.dll.config"
  RMDir "$1\AddIns\OSRAAction"
  no_symyx:
;  strcpy $3 "CambridgeSoft\ChemScript"
;  call un.CheckSoftVersion
;  strcmp $2 "12.0" +1 no_chemoffice
  call un.getChemScriptPath
  strcmp $1 "" no_chemoffice
  Delete "$1\Scripts\Import Structures with OSRA.py"
  no_chemoffice:
SectionEnd

Function CheckSoftVersion
StrCpy $0 0
StrCpy $2 ""
loop:
  EnumRegKey $1 HKLM "Software\$3" $0
  StrCmp $1 "" done
  StrCpy $2 $1
  IntOp $0 $0 + 1
  Goto loop
done:
; $2 contains the version of Soft now or empty
FunctionEnd

Function un.CheckSoftVersion
StrCpy $0 0
StrCpy $2 ""
loop:
  EnumRegKey $1 HKLM "Software\$3" $0
  StrCmp $1 "" done
  StrCpy $2 $1
  IntOp $0 $0 + 1
  Goto loop
done:
; $2 contains the version of Soft now or empty
FunctionEnd

Function downloadPIL
   DetailPrint "need to download and install Python Imaging Library"
   Call ConnectInternet ;Make an internet connection (if no connection available)
   StrCpy $2 "$TEMP\${PIL_EXE}"
   NSISdl::download ${PIL_URL} $2
   Pop $0
   StrCmp $0 success success
    SetDetailsView show
    DetailPrint "download failed: $0"
    Abort
   success:
    ExecWait "$2"
    Delete $2
FunctionEnd

Function getGhostscriptInstPath
 strcmp $2 "" download_gs get_path
 download_gs:
   DetailPrint "need to download and install Ghostscript"
   Call ConnectInternet ;Make an internet connection (if no connection available)
   StrCpy $2 "$TEMP\${GS_EXE}"
   NSISdl::download ${GS_URL} $2
   Pop $0
   StrCmp $0 success success
    SetDetailsView show
    DetailPrint "download failed: $0"
    Abort
   success:
    ExecWait "$2"
    Delete $2
    strcpy $2 "9.07"
	
 get_path:
  strcpy $1 ""
  ReadRegStr $0 HKLM \
     "Software\GPL Ghostscript\$2" \ 
     "GS_DLL"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 -16
  IfFileExists $1\bin\gswin32c.exe fin
  StrCpy $1 ""
  fin:
  ;$1 contains the folder of Ghostscript or empty
FunctionEnd

Function getSymyxPath
 strcpy $1 ""
 ReadRegStr $0 HKLM \
     "Software\Symyx Technologies, Inc.\Symyx Draw\Client\$2" \ 
     "Path"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 
  IfFileExists $1\SymyxDraw.exe fin
  IfFileExists $1\AccelrysDraw.exe fin
  StrCpy $1 ""
  fin:
  ;$1 contains the folder of Symyx Draw or empty
FunctionEnd

Function getChemScriptPath
 strcpy $1 ""
 ReadRegStr $0 HKLM \
     "Software\CambridgeSoft\ChemDraw\12.0\General" \ 
     "ChemDraw Items Default Path"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 -15
  IfFileExists "$1\Scripts\Get 3D Structure.py" fin
  StrCpy $1 ""
  fin:
FunctionEnd

Function un.getSymyxPath
 strcpy $1 ""
 ReadRegStr $0 HKLM \
     "Software\Symyx Technologies, Inc.\Symyx Draw\Client\$2" \ 
     "Path"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 
  IfFileExists $1\SymyxDraw.exe fin
  IfFileExists $1\AccelrysDraw.exe fin
  StrCpy $1 ""
  fin:
  ;$1 contains the folder of Symyx Draw or empty
FunctionEnd

Function un.getChemScriptPath
 strcpy $1 ""
 ReadRegStr $0 HKLM \
     "Software\CambridgeSoft\ChemDraw\12.0\General" \ 
     "ChemDraw Items Default Path"
  StrCmp $0 "" fin extract
  
 extract:
  StrCpy $1 $0 -15
  IfFileExists "$1\Scripts\Get 3D Structure.py" fin
  StrCpy $1 ""
  fin:
FunctionEnd

Function createOSRAbat
fileOpen $0 "$INSTDIR\osra.bat" w
  fileWrite $0 '\
@echo off$\r$\n\
setlocal$\r$\n\
set exec_dir=%~dp0%$\r$\n\
set PATH=%exec_dir%;$1\bin;$1\lib;%PATH%$\r$\n\
set MAGICK_CONFIGURE_PATH=%exec_dir%$\r$\n\
"%exec_dir%osra-bin.exe" %*$\r$\n\
endlocal$\r$\n\
'
fileClose $0
FunctionEnd

Function ReplaceLineStr
 Exch $R0 ; string to replace that whole line with
 Exch
 Exch $R1 ; string that line should start with
 Exch
 Exch 2
 Exch $R2 ; file
 Push $R3 ; file handle
 Push $R4 ; temp file
 Push $R5 ; temp file handle
 Push $R6 ; global
 Push $R7 ; input string length
 Push $R8 ; line string length
 Push $R9 ; global
 
  StrLen $R7 $R1
 
  GetTempFileName $R4
 
  FileOpen $R5 $R4 w
  FileOpen $R3 $R2 r
 
  ReadLoop:
  ClearErrors
   FileRead $R3 $R6
    IfErrors Done
 
   StrLen $R8 $R6
   StrCpy $R9 $R6 $R7 -$R8
   StrCmp $R9 $R1 0 +3
 
    FileWrite $R5 "$R0$\r$\n"
    Goto ReadLoop
 
    FileWrite $R5 $R6
    Goto ReadLoop
 
  Done:
 
  FileClose $R3
  FileClose $R5
 
  SetDetailsPrint none
   Delete $R2
   Rename $R4 $R2
  SetDetailsPrint both
 
 Pop $R9
 Pop $R8
 Pop $R7
 Pop $R6
 Pop $R5
 Pop $R4
 Pop $R3
 Pop $R2
 Pop $R1
 Pop $R0
FunctionEnd


Function ConnectInternet

  Push $R0
    
    ClearErrors
    Dialer::AttemptConnect
    IfErrors noie3
    
    Pop $R0
    StrCmp $R0 "online" connected
      MessageBox MB_OK|MB_ICONSTOP "Cannot connect to the internet."
      Quit
    
    noie3:
  
    ; IE3 not installed
    MessageBox MB_OK|MB_ICONINFORMATION "Please connect to the internet now."
    
    connected:
  
  Pop $R0
  
FunctionEnd


Function .onInit
	# call userInfo plugin to get user info.  The plugin puts the result in the stack
    userInfo::getAccountType
   
    # pop the result from the stack into $0
    pop $0
 
    # compare the result with the string "Admin" to see if the user is admin.
    # If match, jump 3 lines down.
    strCmp $0 "Admin" +3
 
    # if there is not a match, print message and return
    messageBox MB_OK "Please run this with Administrator privileges"
    Quit
 strcpy $3 "Symyx Technologies, Inc.\Symyx Draw\Client"
 call CheckSoftVersion
 strcmp $2 "" no_symyx
 call getSymyxPath
 strcmp $1 "" no_symyx
 SectionGetFlags "${symyx_draw}" $0
 IntOp $0 $0 | ${SF_SELECTED}
 SectionSetFlags "${symyx_draw}" $0
 no_symyx:
 call getChemScriptPath
 strcmp $1 "" no_chemoffice
 SectionGetFlags "${chemoffice}" $0
 IntOp $0 $0 | ${SF_SELECTED}
 SectionSetFlags "${chemoffice}" $0
 no_chemoffice:
FunctionEnd
