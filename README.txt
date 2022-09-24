
To get a brief summary of the important command-line options, invoke lutefisk
with the '-h' option. (On the Macintosh, command-line options are entered on 
the 'Argurments:' line of the initial dialog box. To use command-line options 
under Win32 the executable must be invoked from a DOS prompt such as by using 
the Command Prompt program.)

USAGE:  lutefisk [options] [CID file pathname]

                -o = output file pathname
                -q = quiet mode ON (default OFF)
                -m = precursor ion mass
                -d = details file pathname
                -p = params file pathname
                -r = residues file pathname
                -s = pathnane of file with database sequences to score
                -v = verbose mode ON (default OFF)
                -h = print this help text

________________________________________________________________________

LUTEFISK SOURCE CODE ARCHIVE CONTENTS:
(The Lutefisk_src.tar.gz archive is a tar archive that has been gzipped.)

* Example Lutefisk input file, "QTof_ETYGDMADCCEK.dta"
* Example Lutefisk output file, "QTof_ETYGDMADCCEK.lut"
* This README file, "README"
* The copright file, "COPYRIGHT" 
* The version history, "HISTORY"

* The lutefisk accessory files:
   
    Lutefisk.details
    Lutefisk.params
    Lutefisk.residues
    database.sequence

* C Source code files for Lutefisk:
        Makefile.XXX - Makefiles for compiling Lutefisk for 
                       AIX, IRIX, LINUX, OSF, OSX, or SUN

 LutefiskGlobalDeclarations.c    Declaration of global variables.
 LutefiskMain.c                  Reads parameter, edman, and detail files.
                                  Also contains the function main.
 LutefiskGetAutoTag.c            Finds peptide sequence tags automatically.
 LutefiskGetCID.c                Loads data file and performs peak detection.
 LutefiskMakeGraph.c             Produces N and C terminal sequence graphs.
 LutefiskSummedNode.c            Combines N and C terminal sequence graphs
                                  into a single graph.
 LutefiskHaggis.c		 Finds candidate sequences by searching for
 				  contiguous series of ions that do not necessarily
 				  connect with either N- or C-termini.
 LutefiskSubseqMaker.c           Uses the sequence graph to determine
                                  sequence candidates.
 LutefiskScore.c                 Assigns score and rank to sequence candidates.
 LutefiskXCorr.c                 Performs cross-correlation scoring.
 LutefiskFourier.c               Used for cross-correlation scoring.
 ListRoutines.c                  Handles lists.
 getopt.c                        Interprets line commands when starting the
                                  program (Needed by the Mac and Win32 versions).
 Lutefisk.rsrc                   Resources needed by the Mac version.

 Header files:

 LutefiskDefinitions.h           Contains #defines and struct definitions.
 LutefiskPrototypes.h            Contains function prototypes.
 ListRoutines.h					
 getopt.h			Needed by the Mac and Win32 versions.
________________________________________________________________________

* Compiling on the Macintosh
Current Metrowerks Projects are included in the "Macintosh" folder (provided as a self-
extracting archive to maintain fidelity). If you have an older compiler you will need to 
create a new "Std C Console PPC" project and add the source files as specified above. Set
the preferred heap size to 16M, the minimum heap size to 8M, and the stack size to 512k. You 
may have to change the files to creator 'CWIE' and type 'TEXT' to get Codewarrior to like 
them. Note that in MacOS X the apps can be compiled using the UNIX directions if desired.

* Compiling under Win32
Current Metrowerks Projects are included in the "Win32" folder. If you have an older 
compiler you will need to create a new "C Console App" project and add the source files 
as specified above.

* Compiling on UNIX
After untarring the archive, copy the makefile for your system to "Makefile".
Use the "make lutefisk" command.
________________________________________________________________________

Questions?  Problems?

  Richard S. Johnson
  jsrichar@alum.mit.edu

