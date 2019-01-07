# ReworkAnalysis

Requirements and inclusions:

   THIS LIBARY USES xmltodict.py WHICH WAS PRODUCED BY:

   https://github.com/martinblech/xmltodict.git

   THIS LIBARY IS BUILT ONTOP OF ANACONDA.

   THIS LIBARY ALSO USES pathos FOR MPI AND traitsui for GUI, MUST FIND AND INSTALL FROM:

   https://github.com/uqfoundation/pathos.git

   https://github.com/enthought/traitsui.git

   OR IF USING ANACONDA, JUST USE:

   conda install -c conda-forge pathos=0.2.0

   conda install traitsui

NEW:

  include Setup.py script which uses conda

MAIN PROGRAM IS ./AnalysisGUI.py.

     On first use, a Param.p file is created with some parameters (in pickled format). To initialise this properly follow this steps IN THIS ORDER!

     1. run ./Analysis.py
     2. click on 'Run params' and click Apply on new window.
     3. click on 'General params' and click Apply
     4. exit and reload ./Analysis.py and your good to go.

     Obviously check all the parameters and make sure they are correct!


QUICK TESTS / interactive mode with results:

      running each of the files (with python or ipython if avail):

      python -i ./FlowOpps.py
      python -i ./TwoPtCorrelators.py
      python -i ./ThreePtCorrelators.py
      python -i ./RatioCorrelators.py
      python -i ./FormFactors.py
      python -i ./FFSolve.py

      will run off some default parameter calculation of each command,
      and give you an interactive session with a instance of the class file with all the data contained in it.


COMMON DEBUGGING TIPS:

       Classes are pickled to file. So when they are loaded, they might have older atributes if you have made modifications to the code.
       Avoid this problem by using the optional argument DefWipe=True to force the calculation to occurr.

Notes:

	NOTE on arithmetic operations on classes: Most efficient if lowest level class is on right:

	e.g. SetsOfCorr + Correlator
	     is more efficient than
	     Correlator + SetsOfCorrs


TODO:

	grep TODO *.py to find what needs to be done within the files.

	Params.py Will be reworked. Currently contains default parameters for testing.

	MiscFuncs.py Will hopefully be removed.

	XmlFormatting.py May stay. Taken from old analysis code

	Operator overloading isn't implemented everywhere. Need a clean way of combining self.name

	Make a qmomdir function LatticeParameters (in MomParams.py) to output momentum results in a neat way.
