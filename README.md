# ReworkAnalysis

  This repository uses a custom version of xmltodict.py.

  Run ./Setup.py to download all the conda repositories needed for this package.

Main computing GUI is ./AnalysisGUI.py.

QUICK TESTS / interactive mode with results:

      running each of the files (with python or ipython if avail):

      python -i ./Autocorr.py
      python -i ./BootStrap.py
      python -i ./FlowOpps.py
      python -i ./FlowOpps.py
      python -i ./TwoPtCorrelators.py
      python -i ./ThreePtCorrelators.py
      python -i ./RatioCorrelators.py
      python -i ./FormFactors.py
      python -i ./FFSolve.py

      will run off some default parameter calculation of each command,
      and give you an interactive session with a instance of the class file with all
      the data contained in it.

Plotting manipulation GUI is run by ./ModPlot.py

    Navigate to your pdf file and it will load the pickle file linked to it.
    The final plot will keep the same aspect ratio as the window that shows.

Some other Utilities

    ./CreateTable.py mainly used for reading in plot data and creating a latex or pandas formatted table of data plotted.

    ./ModStylePlot.py uses lists to change the properties of multiple plot planes at once
    (e.g. changing the title size to xyz...)


Notes:

	NOTE on arithmetic operations on classes: Most efficient if lowest level class is on right:

	e.g. SetsOfCorr + Correlator
	     is more efficient than
	     Correlator + SetsOfCorrs

  Just be careful with the saving of functions, they are dumped into a separate folder because
  storing instances of functions to pickle is not recommended

TODO:

	grep TODO *.py to find what needs to be done within the files.

	MiscFuncs.py Will hopefully be removed.

	Operator overloading isn't implemented everywhere. Need a clean way of combining self.name
