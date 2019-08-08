# ReworkAnalysis

see ./build/html/index.html for documentation.

Requirements and inclusions:

   THIS LIBARY USES xmltodict.py WHICH WAS PRODUCED BY:

   https://github.com/martinblech/xmltodict.git

   Run ./Setup.py script which uses conda to setup environment

MAIN PROGRAM IS ./AnalysisGUI.py.

     On first use, a Param.p file is created with some parameters (in pickled format). To initialise this properly follow this steps IN THIS ORDER!

     1. run ./Analysis.py
     2. click on 'Run params' and click Apply on new window.
     3. click on 'General params' and click Apply
     4. exit and reload ./Analysis.py and your good to go.

     Obviously check all the parameters and make sure they are correct!

NEW TEST FILE HAS BEEN ADDED:

      ipython -i ./TestCode.py

      A help message should come up to help you select which test to do.


QUICK TESTS / interactive mode with results:

      running each of the files (with python or ipython if avail):

      ipython -i ./ReadBinaryCfuns.py /path/to/correlator.2cf
      ipython -i ./ReadBinaryCfuns.py /path/to/correlator.3cf
      ipython -i ./ReadFO.py /path/to/flowed_operator.out
      ipython -i ./ReadFO.py /path/to/flowed_operator.out full
      ipython -i ./BootStrapping.py
      ipython -i ./BootStrapping.py
      ipython -i ./Autocorr.py
      ipython -i ./FitFunctions.py
      ipython -i ./FlowOpps.py
      ipython -i ./TwoPtCorrelators.py
      ipython -i ./ThreePtCorrelators.py
      ipython -i ./RatioCorrelators.py
      ipython -i ./FormFactors.py
      ipython -i ./FFSolve.py

      will run off some default parameter calculation of each command,
      and give you an interactive session with a instance of the class file with all the data contained in it.

      Warning! pickling objects in __main__ runs (which is above) will cause issues when unpickling.
      Try to use the actual commands.

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

	MiscFuncs.py Will hopefully be removed.

	XmlFormatting.py May stay. Taken from old analysis code

	Operator overloading isn't implemented everywhere. Need a clean way of combining self.name
