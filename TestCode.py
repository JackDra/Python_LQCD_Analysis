#!/usr/bin/env python

import sys

run_list = [
'BootStrapping',
'Autocorr',
'FitFunctions',
'SetsOfFits',
'FlowOpps',
'TwoPtCorrelators',
'ThreePtCorrelators',
'RatioCorrelators',
'GammaMatricies',
'FormFactors',
'FFSolve'
]

run_list_low = [irun.lower() for irun in run_list]

def HelpPrint():
    print('The allowed input parameters are:')
    for ifile in run_list:
        print('    ',ifile)
    print('Run with command')
    print('ipython -i ./TestCode.py ....')
    print('to obtain the data after the test is complete to test with.')


if len(sys.argv) < 2:
    print('please pass an input parameter to specifiy what file to test')
    HelpPrint()
    sys.exit()

this_run = sys.argv[1].lower()

if this_run not in run_list_low:
    print(sys.argv[1], ' is not in list of files to run')
    HelpPrint()
    sys.exit()


if this_run == 'BootStrapping'.lower():
    from BootStrapping import TestBoot
    data,data2,data3 = TestBoot()
    print('testing data is in data, data2 and data3')
elif this_run == 'Autocorr'.lower():
    from Autocorr import TestAuto
    testdata1,testdata2,testdata3,testdatarat = TestAuto()
    print('Result is in testdata1,2,3')
elif this_run == 'FitFunctions'.lower():
    from FitFunctions import TestComb
    fit_data1,fit_data2,fit_data_comb = TestComb()
    print('data is in fit_data1,fit_data2,fit_data_comb')
elif this_run == 'SetsOfFits'.lower():
    from SetsOfFits import TestSetOfFits,TestSetOfFits_2D
    testdata,testfit = TestSetOfFits()
    testdata_2D,testfit_2D = TestSetOfFits_2D()
    print('data is in testdata,testfit and testdata_2D, testfit_2D')
elif this_run == 'FlowOpps'.lower():
    from FlowOpps import TestFO,TestFOFull
    data,datacomb = TestFO()
    dataFull = TestFOFull()
    print('data is in data,datacomb,dataFull')
elif this_run == 'TwoPtCorrelators'.lower():
    from TwoPtCorrelators import TestTPCorr,TestNNQCorr
    tpdata = TestTPCorr()
    nnqdata = TestNNQCorr()
    print('Testing 2ptcorr complete, see tpdata and nnqdata variables')
elif this_run == 'ThreePtCorrelators'.lower():
    from ThreePtCorrelators import Test3pt,TestNJNQ
    data3 = Test3pt()
    dataNJNQ = TestNJNQ()
    print('Testing 3ptcorr complete, see data3 and dataNJNQ variables')
elif this_run == 'RatioCorrelators'.lower():
    from RatioCorrelators import TestRatFO,TestRat
    dataRatFO = TestRatFO()
    dataRat = TestRat()
    print('Ratio function test complete, see variables dataRat and dataRatFO')
elif this_run == 'FormFactors'.lower():
    from FormFactors import TestFF
    data = TestFF()
    print('data is in data')
elif this_run == 'GammaMatricies'.lower():
    from GammaMatricies import TestSystem
    data = TestSystem()
    print('data is in data')
elif this_run == 'FFSolve'.lower():
    from FFSolve import TestFFSolve
    dataFF = TestFFSolve()
    print('Data is in dataFF')
