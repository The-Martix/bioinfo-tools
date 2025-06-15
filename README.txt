README - How to Import These Libraries from Any Script
=======================================================

To use these libraries (structools, aminotools, seqtools) from any Python script on your system,
you need to add the folder containing them to your global PYTHONPATH. This tells Python where to look for the modules.

-------------------------------------------------------
STEPS BY OPERATING SYSTEM
-------------------------------------------------------

🪟 WINDOWS
----------

1. Copy the absolute path to the folder containing the libraries, for example:
   C:\Users\YourUser\Documents\bioinfo_libs

2. Open a terminal (cmd or PowerShell) and run:

   setx PYTHONPATH "C:\Users\YourUser\Documents\bioinfo_libs"

3. Restart the terminal (or your computer) for the changes to take effect.

🐧 LINUX / 🍎 macOS
------------------

1. Copy the absolute path to the folder, for example:
   /home/your_user/bioinfo_libs

2. Add the following line to the end of your `~/.bashrc`, `~/.zshrc`, or `~/.profile` file:

   export PYTHONPATH="/home/your_user/bioinfo_libs:$PYTHONPATH"

3. Save the file and run:

   source ~/.bashrc       (or the file matching your shell)

-------------------------------------------------------
USAGE IN A SCRIPT
-------------------------------------------------------

Once configured, you can import the modules from anywhere like this:

    from structools import *
    from aminotools import *
    from seqtools import *

-------------------------------------------------------
