from os.path import join
import pickle
from gurobipy import *
from gurobipy import GRB
import scripts.misc_helpers as MH


def ModelGenerator(modeltype, func, extraD):

    with open('../data/common/F-SM-INFO.p', 'rb') as f:
        FIdict = pickle.load(f)
    NQ = FIdict[func]['NQ']
    NK = FIdict[func]['NK']
    NCBS = int(pow(2, NQ))
    MinD = FIdict[func]['MinD']
    ND = MinD + extraD

    CBS = [('{0:0%sb}' % NQ).format(i) for i in range(NCBS)]
    Q = [i for i in range(1, NQ + 1)]
    D = [i for i in range(1, ND + 1)]
    K = [i for i in range(1, NK + 1)]

    with open('../data/common/F-TLPQcoeff.p', 'rb') as f:
        TLPQdict = pickle.load(f)
    with open('../data/pckl/%s.p' % func, 'rb') as f:
        NWRHSdict = pickle.load(f)

    P = TLPQdict[func]
    TARGET_DICT = NWRHSdict['TARGETDICT']
    CBS2COMMODITY_DICT = NWRHSdict['CBS2COMDICT']

    # Model parameter b
    uS, uT, b = {}, {}, {}  # (sigma, k):value of b
    for sigma in CBS:
        uS[sigma] = {}
        for k in K:
            if TARGET_DICT[sigma] == k:
                uS[sigma][k] = 1
            else:
                uS[sigma][k] = 0
    for sigma in CBS:
        uT[sigma] = {}
        for k in K:
            if k in CBS2COMMODITY_DICT[sigma]:
                uT[sigma][k] = 1
            else:
                uT[sigma][k] = 0
    for k in K: b[k] = list(TARGET_DICT.values()).count(k)

    # Set a model
    m = Model()
    x, w, t, xi, Tbar, Wbar = {}, {}, {}, {}, {}, {}

    # Generate variables
    for sigma in CBS:
        for k in K:
            x['S', sigma, k, 0] = m.addVar(vtype=GRB.CONTINUOUS, ub=1.0, name="x_S,%s,K%s,D0" % (sigma, k))
    for d in D:
        for sigma in CBS:
            for k in K:
                x[sigma, sigma, k, d] = m.addVar(vtype=GRB.CONTINUOUS, ub=1.0,
                                                 name="x_%s,%s,K%s,D%s" % (sigma, sigma, k, d))
                for j in Q:
                    pi = MH.Hdist1(sigma, j)
                    x[sigma, pi, k, d] = m.addVar(vtype=GRB.CONTINUOUS, ub=1.0,
                                                  name="x_%s,%s,K%s,D%s" % (sigma, pi, k, d))
    for sigma in CBS:
        for k in K:
            x[sigma, 'T', k, ND + 1] = m.addVar(vtype=GRB.CONTINUOUS, ub=1.0,
                                                name="x_%s,T,K%s,D%s" % (sigma, k, ND + 1))

    for d in D:
        for sigma in CBS:
            xi[sigma, d] = m.addVar(vtype=GRB.BINARY, name="xi_CBS%s,D%s" % (sigma, d))
    for d in D:
        for j in Q:
            w[j, d] = m.addVar(vtype=GRB.BINARY, name="w_Q%s,D%s" % (j, d))
    for d in D:
        for q in range(NQ):
            Wbar[q, d] = m.addVar(vtype=GRB.BINARY, name="Wbar_NC%s,D%s" % (q, d))
    for d in D:
        for j in Q:
            t[j, d] = m.addVar(vtype=GRB.BINARY, name="t_Q%s,D%s" % (j, d))
            t[j, d].BranchPriority = 1
    for j in Q:
        if (P[j]['CHANGE'], P[j]['REMAIN']) != (0, 0):
            for q in [0, 1, 2]:
                Tbar[q, j] = m.addVar(vtype=GRB.BINARY, name="Tbar_NT%s,Q%s" % (q, j))

    CostDict = {3: (1, 1, 5),
                4: (1, 1, 5, 13),
                5: (1, 1, 5, 13, 29),
                6: (1, 1, 5, 13, 29, 61),
                7: (1, 1, 5, 13, 26, 52, 125),
                8: (1, 1, 5, 13, 26, 52, 80, 125)}
    C = quicksum(quicksum(CostDict[NQ][j] * Wbar[j, d] for j in range(NQ)) for d in D)
    m.setObjective(C, GRB.MINIMIZE)
    m.update()

    for k in K:
        m.addConstr(quicksum(x['S', sigma, k, 0] for sigma in CBS) == b[k], name="NWSupply_K%s" % k)
    for sigma in CBS:
        for k in K:
            m.addConstr(
                x['S', sigma, k, 0] == x[sigma, sigma, k, 1] + quicksum(x[sigma, MH.Hdist1(sigma, j), k, 1] for j in Q),
                name="NWConsrv_%s,K%s,L%s" % (sigma, k, 0))
    for d in D[:-1]:
        for sigma in CBS:
            for k in K:
                m.addConstr(x[sigma, sigma, k, d] + quicksum(x[MH.Hdist1(sigma, j), sigma, k, d] for j in Q) - x[sigma, sigma, k, d + 1]
                            - quicksum(x[sigma, MH.Hdist1(sigma, j), k, d + 1] for j in Q) == 0,
                            name="NWConsrv_%s,K%s,L%s" % (sigma, k, d))
    for sigma in CBS:
        for k in K:
            m.addConstr(x[sigma, sigma, k, ND] + quicksum(x[MH.Hdist1(sigma, j), sigma, k, ND] for j in Q) == x[
                sigma, 'T', k, ND + 1],
                        name="NWConsrv_%s,K%s,L%s" % (sigma, k, ND))
    for k in K:
        m.addConstr(quicksum(x[sigma, 'T', k, ND + 1] for sigma in CBS)
                    == b[k],
                    name="NWDemand_K%s" % (k))

    for sigma in CBS:
        for k in K:
            m.addConstr(x['S', sigma, k, 0]
                        == uS[sigma][k],
                        name="NW,UB_S,%s,K%s,D%s" % (sigma, k, 0))

    for sigma in CBS:
        for k in K:
            if uT[sigma][k] == 0:
                m.addConstr(x[sigma, 'T', k, ND + 1] <= uT[sigma][k], name="NW,UB_%s,T,K%s,D%s" % (sigma, k, ND + 1))
    for d in D:
        for sigma in CBS:
            m.addConstr(quicksum(x[sigma, sigma, k, d] for k in K) == xi[sigma, d],
                        name="NWCIR1_%s,%s,D%s" % (sigma, sigma, d))
            for j in Q:
                pi = MH.Hdist1(sigma, j)
                m.addConstr(quicksum(x[sigma, pi, k, d] for k in K) <= 1 - xi[sigma, d],
                            name="NWCIR2a_%s,%s,Q%s,D%s" % (sigma, pi, j, d))
                m.addConstr(quicksum(x[sigma, pi, k, d] for k in K) <= t[j, d],
                            name="NWCIR2b_%s,%s,Q%s,D%s" % (sigma, pi, j, d))
                m.addConstr(quicksum(x[sigma, pi, k, d] for k in K) >= t[j, d] - xi[sigma, d],
                            name="NWCIR2c_%s,%s,Q%s,D%s" % (sigma, pi, j, d))

    for d in D:
        for sigma in CBS:
            SigmaDict = MH.SigmaDict(sigma, NQ)
            tmpM = len([j for j in Q if SigmaDict[j] == 0])
            m.addConstr(xi[sigma, d] <= quicksum(w[j, d] for j in Q) - quicksum(w[j, d] * SigmaDict[j] for j in Q),
                        name="CFTEST,L_CBS%s,D%s" % (sigma, d))
            try:
                if modeltype == 'MS':
                    for j in Q:
                        if SigmaDict[j] == 0:
                            m.addConstr(w[j, d] <= xi[sigma, d], name="CFTEST,R_CBS%s,D%s" % (sigma, d))
                elif modeltype == 'MW':
                    m.addConstr(quicksum(w[j, d] for j in Q) - quicksum(w[j, d] * SigmaDict[j] for j in Q) <= tmpM * xi[sigma, d], name="CFTEST,R_CBS%s,D%s" % (sigma, d))
            except:
                raise ("Model type error: None of MS or MW")



    for d in D:
        for j in Q:
            m.addConstr(t[j, d] + w[j, d] <= 1, name="CIR,TW_Q%s,D%s" % (j, d))
    for d in D:
        m.addConstr(quicksum(t[j, d] for j in Q) == 1, name="CIR,T_D%s" % (d))
    for d in D:
        m.addConstr(quicksum(j * Wbar[j, d] for j in range(NQ)) - quicksum(w[j, d] for j in Q) == 0, name="QCa_D%s" % (d))
        m.addConstr(quicksum(Wbar[j, d] for j in range(NQ)) == 1, name="QCb_D%s" % (d))
    for j in Q:
        if (P[j]['CHANGE'], P[j]['REMAIN']) != (0, 0):
            m.addConstr(quicksum(t[j, d] for d in D) <= Tbar[1, j] + ND * Tbar[2, j])
            m.addConstr(Tbar[1, j] + 2 * Tbar[2, j] <= quicksum(t[j, d] for d in D))
            m.addConstr(quicksum(Tbar[q, j] for q in [0, 1, 2]) == 1)
            m.addConstr(Tbar[0, j] <= 1 - P[j]['CHANGE'] + P[j]['REMAIN'])
            m.addConstr(Tbar[1, j] <= P[j]['CHANGE'] + 1 - P[j]['REMAIN'])

    try:
        PRJpath = "../models"
        MODdir = join(PRJpath, "%s-%s-D%s+%s-test.mps" % (modeltype, func, MinD, extraD))
        m.write(MODdir)

    except Exception as e:
        print("Error reported while saving the model as an mps file: ", e)