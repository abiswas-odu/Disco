BBMap/BBTools readme
Written by Brian Bushnell
Last updated December 23, 2015

The BBTools package was written by Brian Bushnell, with the exception of the (optional, but faster) C, JNI, and MPI components, which were written by Jonathan Rood.

All tools in the BBTools package are free to use.  If you use BBTools in work leading to a publication, and BBTools has not yet been published, please cite it something like this:
BBMap - Bushnell B. - sourceforge.net/projects/bbmap/

License:

The BBMap package is open source and free to use with no restrictions.  For more information, please read Legal.txt and license.txt.

Documentation:

Documentation is in the /bbmap/docs/ directory, and in each tool's shellscript in /bbmap/.
readme.txt: This file.
UsageGuide.txt: Contains basic installation and usage information.  Please read this first!
ToolDescriptions.txt: Contains a list of all BBTools, a description of what they do, and their hardware requirements.
compiling.txt: Information on compiling JNI code.
readme_config.txt: Usage information about config files.
readme_filetypes.txt: More detailed information on file formats supported by BBTools.
changelog.txt: List of changes by version, and current known issues.

Tool-specific Guides:

Some tools have specific guides, like BBDukGuide.txt.  They are in /bbmap/docs/guides/.  For complete documentation of a tool, I recommend that you read UsageGuide.txt first (which covers the shared functionality of all tools), then the tool's specific guide if it has one (such as ReformatGuide.txt), then the tool's shellscript (such as reformat.sh) which lists all of the flags.

If you have any questions not answered in the documentation, please look at the relevant SeqAnswers thread (linked from here: http://seqanswers.com/forums/showthread.php?t=41057) and post a question there if it is not already answered.  You can also contact JGI's BBTools team at bbtools@lbl.gov, or me at bbushnell@lbl.gov.  But please read the documentation first.

Special thanks for help with shellscripts goes to:
Alex Copeland (JGI), Douglas Jacobsen (JGI/NERSC), Bill Andreopoulos (JGI), sdriscoll (SeqAnswers), Jon Rood (JGI/NERSC), and Elmar Pruesse (UC Denver).

Special thanks for helping to support BBTools goes to Genomax (SeqAnswers).
