import os
import pickle
import argparse
from gurobipy import *
from utils import *
func, extra_d, env_id = None, None, None


def exp_setting():
    global func, extra_d, env_id
    parser = argparse.ArgumentParser(description="Set func, extra_d, and env_id")

    # Adding arguments
    parser.add_argument('--func', type=str, help="Function name")
    parser.add_argument('--extra_d', type=int, help="Extra gates")
    parser.add_argument('--env_id', type=str, help="Configuration ID")

    # Parsing arguments
    args = parser.parse_args()

    # Accessing variables
    func = args.func
    extra_d = args.extra_d
    env_id = args.env_id

    # Print or use the variables as needed
    print(f"Function: {func}")
    print(f"Extra Gates: {extra_d}")
    print(f"Environment ID: {env_id}")


if __name__ == '__main__':

    # *****************
    # EXP SETTING
    # *****************
    exp_setting()

    # *****************
    # SOLVER SETTING
    # *****************
    # key: evn_id (str)(E1~E6/E1T~E6T)
    # value = tuple of length 4 or 5 (Optional)
    # value[0] -> Model Type (str)(MW/MS)
    # value[1] -> Branching Priority (Bool)
    # value[2] -> GRB cuts Param (int)(0/3)
    # value[3] -> USER xi-cuts Param (Bool)
    # value[4] -> BnC Alg. Key (str)(BnC/CnB)
    env_dict = {'E1': ('MW', False, 0, False),
                'E2': ('MW', False, 3, False),
                'E3': ('MW', False, 0, True, 'BnC'),
                'E4': ('MW', False, 0, True, 'CnB'),
                'E5': ('MS', False, 0, False),
                'E6': ('MS', False, 3, False),
                'E1T': ('MW', True, 0, False),
                'E2T': ('MW', True, 3, False),
                'E3T': ('MW', True, 0, True, 'BnC'),
                'E4T': ('MW', True, 0, True, 'CnB'),
                'E5T': ('MS', True, 0, False),
                'E6T': ('MS', True, 3, False),
                }

    default_param = [("TimeLimit", 7200),
                     ("Method", 2),
                     ("NodeMethod", 1),
                     ("DegenMoves", 0),
                     ("SubMIPNodes", 20),
                     ("Presolve", 0),
                     ("Precrush", 1),
                     ("Heuristics", 0), ]

    model_type, bp_key, grb_param_cut, xicut_key = (env_dict[env_id][0],
                                                    env_dict[env_id][1],
                                                    env_dict[env_id][2],
                                                    env_dict[env_id][3])


    # *****************
    # DATA RETRIEVAL
    # *****************
    curr_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(os.getcwd(), 'data/common/func_attrs.p'), 'rb') as f:
        FIdict = pickle.load(f)
    min_d = FIdict[func]['MinD']

    NQ = FIdict[func]['NQ']
    NK = FIdict[func]['NK']
    ND = min_d + extra_d
    NCBS = int(pow(2, NQ))

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


    # *****************
    # MODEL RECOVERY
    # *****************

    m = read(os.path.join(os.getcwd(), 'models', '%s-%s-D%s+%s.mps' % (model_type, func, min_d, extra_d)))
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
    for k in K:
        b[k] = list(target_dict.values()).count(k)

    x, w, t, xi, Tbar, Wbar = {}, {}, {}, {}, {}, {}
    for sigma in CBS:
        for k in K:
            x['S', sigma, k, 0] = m.getVarByName("x_S,%s,K%s,D0" % (sigma, k))
    for d in D:
        for sigma in CBS:
            for k in K:
                x[sigma, sigma, k, d] = m.getVarByName("x_%s,%s,K%s,D%s" % (sigma, sigma, k, d))
                for j in Q:
                    pi = add_unit_hdist(sigma, j)
                    x[sigma, pi, k, d] = m.getVarByName("x_%s,%s,K%s,D%s" % (sigma, pi, k, d))
    for sigma in CBS:
        for k in K:
            x[sigma, 'T', k, ND + 1] = m.getVarByName("x_%s,T,K%s,D%s" % (sigma, k, ND + 1))
    for d in D:
        for sigma in CBS:
            xi[sigma, d] = m.getVarByName("xi_CBS%s,D%s" % (sigma, d))
    for d in D:
        for j in Q:
            w[j, d] = m.getVarByName("w_Q%s,D%s" % (j, d))
    for d in D:
        for q in range(NQ):
            Wbar[q, d] = m.getVarByName("Wbar_NC%s,D%s" % (q, d))
    for d in D:
        for j in Q:
            t[j, d] = m.getVarByName("t_Q%s,D%s" % (j, d))
    for j in Q:
        if (P[j]['CHANGE'], P[j]['REMAIN']) != (0, 0):
            for q in [0, 1, 2]:
                Tbar[q, j] = m.getVarByName("Tbar_NT%s,Q%s" % (q, j))


    # *****************
    # CALLBACKS
    # *****************

    def callback_default(model, where):
        if where == GRB.Callback.MIP:
            crt_cutcnt = model.cbGet(GRB.Callback.MIP_CUTCNT)
            if crt_cutcnt != model._cutcnt:
                print('*** %s cut(s) added...' % crt_cutcnt)
                model._cutcnt = crt_cutcnt


    def callback_cnb(model, where):
        global D
        global CBS

        if where == GRB.Callback.MIP:
            crt_cutcnt = model.cbGet(GRB.Callback.MIP_CUTCNT)
            if crt_cutcnt != model._cutcnt:
                print('*** %s cut(s) added...' % crt_cutcnt)
                model._cutcnt = crt_cutcnt
        if where == GRB.Callback.MIPNODE:
            if model.cbGet(GRB.Callback.MIPNODE_NODCNT) == 0:
                sol_w = model.cbGetNodeRel(model._w)
                sol_xi = model.cbGetNodeRel(model._xi)
                for d in D:
                    for sigma in CBS:
                        temp_sigma = create_sigma_dict(sigma, NQ)
                        for j in Q:
                            if temp_sigma[j] == 0:
                                if sol_xi[sigma, d] <= sol_w[j, d]:
                                    model.cbCut(model._xi[sigma, d] >= model._w[j, d])
                                    model._xicut += 1


    def callback_bnc(model, where):
        global D
        global CBS

        if where == GRB.Callback.MIP:
            crt_cutcnt = model.cbGet(GRB.Callback.MIP_CUTCNT)
            if crt_cutcnt != model._cutcnt:
                print('*** %s cut(s) added...' % crt_cutcnt)
                model._cutcnt = crt_cutcnt
        if where == GRB.Callback.MIPNODE:
            if model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
                sol_w = model.cbGetNodeRel(model._w)
                sol_xi = model.cbGetNodeRel(model._xi)
                for d in D:
                    for sigma in CBS:
                        temp_sigma = create_sigma_dict(sigma, NQ)
                        for j in Q:
                            if temp_sigma[j] == 0:
                                if sol_xi[sigma, d] <= sol_w[j, d]:
                                    model.cbCut(model._xi[sigma, d] >= model._w[j, d])
                                    model._xicut += 1


    # *****************
    # SOLVING & LOGGING
    # *****************

    log_path = os.path.join(os.getcwd(), "results/log", "%s-D%s+%s-%s-TEST.txt" % (func, min_d, extra_d, env_id))
    sol_path = os.path.join(os.getcwd(), "results/sol", "%s-D%s+%s-%s-TEST.sol" % (func, min_d, extra_d, env_id))

    log_start(log_path)
    for param in default_param:
        m.setParam(param[0], param[1])
    m.setParam("Cuts", grb_param_cut)

    if xicut_key:
        xicut_app = env_dict[env_id][4]
    else:
        xicut_app = None

    if bp_key:
        for key in list(t.keys()):
            t[key[0], key[1]].BranchPriority = 1

    m._cutcnt = 0
    m._xicut = 0
    if xicut_key:
        m._w = w
        m._xi = xi
        if xicut_app == 'BnC':
            m.optimize(callback_bnc)
        elif xicut_app == 'CnB':
            m.optimize(callback_cnb)
    else:
        m.optimize(callback_default)
    log_stop()

    sol_dict = create_final_sol_dict(m, x, t, w, xi, Wbar, Tbar)
    with open(sol_path, 'wb') as f:
        pickle.dump(sol_dict, f)