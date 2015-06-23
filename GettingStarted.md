# Running Polymode on Windows #

## Installing ##

**Polymode** requires the **Python(x,y)** scientific Python environment to be installed. It is available from
http://www.pythonxy.com/

First download the **Python(x,y)** web installer available on the  [Python(x,y) website](http://pythonxy.com/dl.php?file=windows/xyweb-1.0.3.exe), execute the installer to start.

Select the light version and _download_. The packages will be downloaded to your computer and the installation will begin.

![http://polymode.googlecode.com/svn/wiki/figures/gs_selectversion.png](http://polymode.googlecode.com/svn/wiki/figures/gs_selectversion.png)

When the installer starts it will ask you to choose optional components. Make sure that you select either **All users** or **Just me** under the "Install for" option. If you select **Just me** then only you will be able to use Python(X,Y).

![http://polymode.googlecode.com/svn/wiki/figures/gs_installfor.png](http://polymode.googlecode.com/svn/wiki/figures/gs_installfor.png)

After the installation is finished Polymode can be installed by running the Polymode installer.
[Polymode website](http://code.google.com/p/polymode/downloads/list)

![http://polymode.googlecode.com/svn/wiki/figures/gs_installfor.png](http://polymode.googlecode.com/svn/wiki/figures/gs_installfor.png)

You now have everything you need to start.

## Using the Python(x,y) Environment ##

The **IPython** command prompt is used to run all code.
This is installed with Python(x,y), you can access it through the Windows Start Menu under
Python (x,y) -> Enhanced Consoles -> IPython(x,y) as shown below:

![http://polymode.googlecode.com/svn/wiki/figures/gs_pythonxy.png](http://polymode.googlecode.com/svn/wiki/figures/gs_pythonxy.png)

You can drag this icon to the main start menu to access it faster:

![http://polymode.googlecode.com/svn/wiki/figures/gs_pythonxy_quick.png](http://polymode.googlecode.com/svn/wiki/figures/gs_pythonxy_quick.png)

<a href='Hidden comment: 
Both an editor *SciTE* and the *IPython* shell can be started using this command window, under the applications menu different applications can be selected and run. The most useful is *SciTE* which allows you to edit python files easily. You can run it selecting the tick button next to the menu.
From the "Interactive consoles" select *IPython(x,y)* and click the leftmost console button, as shown. This will start the IPython shell in a new window.
'></a>

In the IPython window you can change directories to where your code is using the `cd` command, `ls` will give you a list of files in the current directory.

The **tab** key can be used to auto-complete any partially completed commands. This can be helpful going through directories and running commands.

A python script file can be run by using the run command, to run a file named "ex1\_materials.py" just enter:

```
In [1]: run ex1_materials
```

Also any command's help can be accessed by typing the command then a \u201c?\u201d afterwards, for example:

```
In [2]: sin? 
Type:           ufunc 
Base Class:     <type 'numpy.ufunc'> 
String Form:    <ufunc 'sin'> 
Namespace:      Interactive 
Docstring: 
    y = sin(x) sine elementwise. 
```

Some other useful commands are:

  * `ls`: List the files in the current directory
  * `pwd`: Return the current working directory path.
  * `edit`: starts an editor and allows you to edit the file. After closing the window the file will be run.
  * `history`: Print input history with most recent last.
  * `logstart`: Log an interactive session to a file
  * `quickref`: Show a quick reference sheet
  * `rep`: Repeat previous input.
  * `store`: Quickly store variables in a file
  * `whos`: List used variables in the workspace
  * `bookmark`: Bookmark a folder
  * `logon`/`logoff`: Logs interactive commands to a file

Help for all of these commands is available in IPython by typing the command then a "?" and hitting enter

### Registry Entries in MPOF Computers ###

The University of Sydney Physics (MPOF) computers do not respect the local system registry entries, to correct this download this registry file and double click on it. Click "	yes" when asked if you want to add the information to the registry.

[Download Python(x,y) registry file](http://polymode.googlecode.com/svn/misc/mpof_pythonxy_registry.reg)

### Paths in Windows ###

The **PATH** environment variable needs to be set correctly for Python(x,y) to run properly and find the required
programs. This is automatically set up when Python(x,y) is installed, however if another user wants to use the
same installation then each user needs to set up the environment variables

In particular the **PATH** variable should include these directories:
```
C:\Python25;C:\Python25\DLLs;C:\Python25\Scripts;C:\Python25\Lib\site-packages\PyQt4;C:\PythonXY\console;C:\PythonXY\SciTE-1.77.1
```

See the [Microsoft support page](http://support.microsoft.com/kb/310519) for information on how to change the path.

### Advanced: Compiling on Windows ###

To compile the code (from svn for example) on Windows machines you need to have quite a few additional dependencies. The bare essentials are Python(x,y) as explained above and a compiler, such as MinGW.

To install [MinGW](http://www.mingw.org) follow the instructions here in the
[Howto install MinGW](http://www.mingw.org/wiki/HOWTO_Install_the_MinGW_GCC_Compiler_Suite)
page on their website. Install compilers for c++ and fortran. Also needed are the ATLAS binaries, either
download them or build them according to the [Numpy/Scipy how to](http://www.scipy.org/Installing_SciPy/Windows). Use the Numpy `site.cfg` file in the `other` directory and customize it to point to the ATLAS libraries as detailed in the how to.

Finally, build Polymode using the command:
```
python.exe setup.py config --compiler=mingw32 build --compiler=mingw32 bdist_wininst
```

_This currently doesn't work as the block lu solver requires Boost, it will be fixed soon hopefully!_

## Installing on Linux/Mac OSX ##

For Ubuntu Jaunty Jackalope (9.04) and later with Python 2.6 installed there are binary downloads available on the `downloads` page.

Some prerequisites are required before Polymode can be installed. These are the python interpreter, an interactive python shell (IPython), numerical python (numpy and scipy) and the Matplotlib plotting library (matplotlib). These can be installed automatically on most linux distributions.

On Ubuntu execute the command:
```
sudo aptitude install ipython python-matplotlib python-numpy python-scipy python-setuptools
```

To enable the block LU solver on Ubuntu Karmic (9.10) install the boost libraries:
```
sudo aptitude install libboost-python1.38-0
```

Finally, download the appropriate egg file for your distribution (32 bit or 64 bit) and install with Setuptools:
```
sudo easy_install [name of egg file]
```

### Building from source ###

To compile the code from source, you'll need to install some development packages. In particular a fortran compiler, Atlas libraries, numpy and Python development packages and boost c++ development libraries.

On Ubuntu these requirements can be installed with the command:
```
sudo aptitude install build-essential gfortran libatlas-sse-dev subversion python-dev
                      libboost-dev
```

Get the software source from the repository:
```
cd <path_to_source_dir>
svn checkout http://polymode.googlecode.com/svn/trunk/ polymode
```

And now you can compile the package using:
```
cd <path_to_source_dir>/polymode
python setup.py build
sudo python setup.py install
```

This should build and install the package for you, if the prerequisites are also installed.

Alternatively an Setuptools egg package can be built with the command:
```
python setup.py bdist_egg
```

### Installing on Mac OS X ###

Currently there is no compiled package for OS X. Numpy, scipy and matplotlib are available in the form of an all-in one installer on the website http://trichech.us/?page_id=5.

At this time the package needs to be compiled from source.

### Installing on the APAC Altix Cluster computers ###
Precompiled versions of the prerequisites for running computing jobs on the APAC Alitx cluster computers are available. Copy these files and the install file to the APAC computers and run the "install.sh" file. It will ask you for a directory and install the required files to it. You will need to add the some commands to your profile file to complete the installation, the installer will tell you what to add.