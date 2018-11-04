import FFSolve as ffs


def ReadFFTop(thisInfo):

    data = ffs.FFEquation('VectorTop',Info=thisInfo)
    data.LoadPickle(DefWipe=False,HardWipe=False)
