import pickle
from os.path import join
from gurobipy import *

import sys
class Transcript(object):

    def __init__(self, filename):
        self.terminal = sys.stdout
        self.logfile = open(filename, "a")

    def write(self, message):
        self.logfile.write(message)

    def flush(self):
        pass

def log_start(filename):
    """Start transcript, appending print output to given filename"""
    sys.stdout = Transcript(filename)

def log_stop():
    """Stop transcript and return print functionality to normal"""
    sys.stdout.logfile.close()
    sys.stdout = sys.stdout.terminal


def Hdist1(single_CBS, j):
    signle_CBS_list = list(single_CBS)
    if signle_CBS_list[j-1] == '0':
        signle_CBS_list[j-1] = '1'
        output_CBS = ''.join(signle_CBS_list)
    else:
        signle_CBS_list[j-1] = '0'
        output_CBS = ''.join(signle_CBS_list)
    return output_CBS


def SigmaDict(sigma, Q_num):
    index_list = [i+1 for i in range(Q_num)]
    bit_list = [int(s) for s in list(sigma)]
    sigma_dict = dict(zip(index_list, bit_list))
    return sigma_dict


def TmpSolDict(GRBdict):
    SolDict = {}
    keylist = list(GRBdict.keys())
    for key in keylist:
        SolDict[key] = GRBdict[key].x
    return SolDict

### [S] Integrate TmpSoldict to the FinalSolDict
def FinalSolDict(m, x, t, w, xi, Wbar, Tbar):
    VarClass = ['Z', 'x', 't', 'w', 'xi', 'Wbar', 'Tbar']
    SolDict = {}
    for V in VarClass:
        if V == 'x': SolDict[V] = TmpSolDict(x)
        if V == 't': SolDict[V] = TmpSolDict(t)
        if V == 'w': SolDict[V] = TmpSolDict(w)
        if V == 'xi': SolDict[V] = TmpSolDict(xi)
        if V == 'Wbar': SolDict[V] = TmpSolDict(Wbar)
        if V == 'Tbar': SolDict[V] = TmpSolDict(Tbar)
        if V == 'Z': SolDict[V] = m.ObjVal
    return SolDict





def decimal2binary(Ndec, Q):
    Nbin = bin(Ndec)[2:].zfill(Q)
    return Nbin


def generateInitState(Q):
    initState = {}
    for i in range(int(pow(2, Q))):
        initState[i] = decimal2binary(i, Q)
    return initState


def bitSwitch(inputCBS, switchIndex):
    outputCBS = list(inputCBS)
    if inputCBS[switchIndex] == '1':
        outputCBS[switchIndex] = '0'
    else:
        outputCBS[switchIndex] = '1'
    return ("".join(outputCBS))


def flipTarget(gate, inputCBS, Q):
    for i in range(Q):
        if gate[i] == 'T':
            outputCBS = bitSwitch(inputCBS, i)
    return outputCBS

def controlSyncTest(gate, inputCBS ,Q):
    Cnum = gate.count('C')
    tmpCNum = 0
    for i in range(Q):
        if gate[i] == 'C' and inputCBS[i] == '1':
            tmpCNum += 1
    if tmpCNum == Cnum:
        return True
    else:
        return False


def generateOutputCBSdict(gate, inputState, Q):
    outputCBSdict = {}
    if gate == 'E'*Q:
        for c in range(pow(2, Q)):
            outputCBSdict[c] = inputState[c]
    else:
        for c in range(pow(2, Q)):
            if controlSyncTest(gate, inputState[c], Q) == True:
                outputCBSdict[c] = flipTarget(gate, inputState[c], Q)
            else:
                outputCBSdict[c] = inputState[c]
    return outputCBSdict


def compare2truthtable(outputCBS, truthtable, Q):
    TFresult, errorBit = True, []
    for i in range(Q):
        if truthtable[i] == '-':
            continue
        elif truthtable[i] != outputCBS[i]:
            TFresult = False
            errorBit.append(i)
            continue
        else:
            continue
    return TFresult, errorBit


def funcSelect(NQrange, NDrange, IC):
    with open('F-SM-INFO.p', 'rb') as f:
        InfoDict = pickle.load(f)

    Qfunc = []
    for key in list(InfoDict.keys()):
        if InfoDict[key]['NQ'] in NQrange:
            Qfunc.append(key)
    ICfunc = []
    for key in list(InfoDict.keys()):
        if InfoDict[key]['IC'] in IC:
            ICfunc.append(key)
    MinDfunc = []
    for key in list(InfoDict.keys()):
        if InfoDict[key]['MinD'] in NDrange:
            MinDfunc.append(key)

    funcList = list(set(Qfunc) & set(ICfunc) & set(MinDfunc))
    return funcList


def funcSizeDict():
    FsizeDict = {'S':[], 'M':[], 'L':[]}
    with open('F-SM-INFO.p', 'rb') as f:
        InfoDict = pickle.load(f)
    for key in list(InfoDict.keys()):
        NQ = InfoDict[key]['NQ']
        MinD = InfoDict[key]['MinD']
        if NQ*MinD <= 24:
            FsizeDict['S'].append(key)
        elif 15 <= NQ*MinD <= 30:
            FsizeDict['M'].append(key)
        else:
            FsizeDict['L'].append(key)
    return FsizeDict


def generateEXPset(Flist, Dlist, Slist):
    exps = []
    for F in Flist:
        for D in Dlist:
            for S in Slist:
                exps.append((F, D, S))
    expind = list(range(1, 1 + len(exps)))
    EXPset = dict(zip(expind, exps))
    return EXPset


def generateEXPfolders(EXPdir):
    nested = ['LOG', 'PLOG', 'SOL', 'CHK', 'CIR']
    nesteddirs = [join(EXPdir, f) for f in nested]
    if not os.path.exists(EXPdir):
        os.makedirs(EXPdir)
        for dir in nesteddirs:
            os.makedirs(dir)


def printExpInfo(EXPsets, key):
    with open('F-SM-INFO.p', 'rb') as f:
        InfoDict = pickle.load(f)
    F, D, S = EXPsets[key][0], EXPsets[key][1], EXPsets[key][2]
    MinD = InfoDict[F]['MinD']
    NEXP = len(list(EXPsets.keys()))
    print("*************************")
    print("(%s/%s) %s-D%s+%s-%s" % (key, NEXP, F, MinD, D, S) )
    print("    1. NQ: %s" % (InfoDict[F]['NQ']))
    print("    2. MinD: %s" % (InfoDict[F]['MinD']))
    print("    3. IC: %s" % (InfoDict[F]['IC']))
    print("*************************")


def solDict2gateDict(solDict, D, Q):
    gateDict = {}
    for i in range(1, D + 1):
        tmpdict = {}
        for j in range(1, Q + 1):
            tmpdict[j] = ''
        gateDict[i] = tmpdict

    for i in range(1, D + 1):
        for j in range(1, Q + 1):
            if (i, j) in list(solDict.keys()):
                gateDict[i][j] += solDict[(i, j)]
            else:
                gateDict[i][j] += 'E'
    return gateDict