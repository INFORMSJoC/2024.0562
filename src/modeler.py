import os
import pickle
import argparse
from utils import *
from gurobipy import *
_model_type, _func, _extra_d = None, None, None


def model_setting():
    global _model_type, _func, _extra_d
    parser = argparse.ArgumentParser(description="Set model type, function, and extra gates")

    # Adding arguments
    parser.add_argument('--func', type=str, help="Function name")
    parser.add_argument('--model_type', type=str, help="Extra gates")
    parser.add_argument('--extra_d', type=int, help="Configuration ID")

    # Parsing arguments
    args = parser.parse_args()

    # Accessing variables
    _model_type = args.model_type
    _func = args.func
    _extra_d = args.extra_d

    # Print or use the variables as needed
    print(f"Model type: {_model_type}")
    print(f"Function: {_func}")
    print(f"Extra Gates: {_extra_d}")


def model_generator(model_type: str, func: str, extra_d: int):
    with open(os.path.join(os.getcwd(), 'data/common/func_attrs.p'), 'rb') as f:
        FIdict = pickle.load(f)
    NQ = FIdict[func]['NQ']
    NK = FIdict[func]['NK']
    NCBS = int(pow(2, NQ))
    min_d = FIdict[func]['MinD']
    ND = min_d + extra_d

    CBS = [('{0:0%sb}' % NQ).format(i) for i in range(NCBS)]
    Q = [i for i in range(1, NQ + 1)]
    D = [i for i in range(1, ND + 1)]
    K = [i for i in range(1, NK + 1)]

    with open(os.path.join(os.getcwd(), 'data/common/func_coeffs.p'), 'rb') as f:
        psp_coeffs_dict = pickle.load(f)
    with open(os.path.join(os.getcwd(), 'data/processed/%s.p') % func, 'rb') as f:
        proc_data_dict = pickle.load(f)

    P = psp_coeffs_dict[func]
    target_dict = proc_data_dict['TARGETDICT']
    cbs2com_dict = proc_data_dict['CBS2COMDICT']

    # Model parameter b
    uS, uT, b = {}, {}, {}  # (sigma, k):value of b
    for sigma in CBS:
        uS[sigma] = {}
        for k in K:
            if target_dict[sigma] == k:
                uS[sigma][k] = 1
            else:
                uS[sigma][k] = 0
    for sigma in CBS:
        uT[sigma] = {}
        for k in K:
            if k in cbs2com_dict[sigma]:
                uT[sigma][k] = 1
            else:
                uT[sigma][k] = 0
    for k in K: b[k] = list(target_dict.values()).count(k)

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
                    pi = add_unit_hdist(sigma, j)
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
                x['S', sigma, k, 0] == x[sigma, sigma, k, 1] + quicksum(x[sigma, add_unit_hdist(sigma, j), k, 1] for j in Q),
                name="NWConsrv_%s,K%s,L%s" % (sigma, k, 0))
    for d in D[:-1]:
        for sigma in CBS:
            for k in K:
                m.addConstr(x[sigma, sigma, k, d] + quicksum(x[add_unit_hdist(sigma, j), sigma, k, d] for j in Q) - x[sigma, sigma, k, d + 1]
                            - quicksum(x[sigma, add_unit_hdist(sigma, j), k, d + 1] for j in Q) == 0,
                            name="NWConsrv_%s,K%s,L%s" % (sigma, k, d))
    for sigma in CBS:
        for k in K:
            m.addConstr(x[sigma, sigma, k, ND] + quicksum(x[add_unit_hdist(sigma, j), sigma, k, ND] for j in Q) == x[
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
                pi = add_unit_hdist(sigma, j)
                m.addConstr(quicksum(x[sigma, pi, k, d] for k in K) <= 1 - xi[sigma, d],
                            name="NWCIR2a_%s,%s,Q%s,D%s" % (sigma, pi, j, d))
                m.addConstr(quicksum(x[sigma, pi, k, d] for k in K) <= t[j, d],
                            name="NWCIR2b_%s,%s,Q%s,D%s" % (sigma, pi, j, d))
                m.addConstr(quicksum(x[sigma, pi, k, d] for k in K) >= t[j, d] - xi[sigma, d],
                            name="NWCIR2c_%s,%s,Q%s,D%s" % (sigma, pi, j, d))

    for d in D:
        for sigma in CBS:
            sigma_dict = create_sigma_dict(sigma, NQ)
            tmpM = len([j for j in Q if sigma_dict[j] == 0])
            m.addConstr(xi[sigma, d] <= quicksum(w[j, d] for j in Q) - quicksum(w[j, d] * sigma_dict[j] for j in Q),
                        name="CFTEST,L_CBS%s,D%s" % (sigma, d))
            try:
                if model_type == 'MS':
                    for j in Q:
                        if sigma_dict[j] == 0:
                            m.addConstr(w[j, d] <= xi[sigma, d], name="CFTEST,R_CBS%s,%s,D%s" % (sigma, j, d))
                elif model_type == 'MW':
                    m.addConstr(quicksum(w[j, d] for j in Q) - quicksum(w[j, d] * sigma_dict[j] for j in Q) <= tmpM * xi[sigma, d], name="CFTEST,R_CBS%s,D%s" % (sigma, d))
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
        MODdir = os.path.join(os.getcwd(), "models", "%s-%s-D%s+%s-test.mps" % (model_type, func, min_d, extra_d))
        m.write(MODdir)

    except Exception as e:
        print("Error reported while saving the model as an mps file: ", e)


if __name__ == '__main__':
    model_setting()
    model_generator(_model_type, _func, _extra_d)
