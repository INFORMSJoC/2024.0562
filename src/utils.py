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
    sys.stdout = Transcript(filename)


def log_stop():
    sys.stdout.logfile.close()
    sys.stdout = sys.stdout.terminal


def add_unit_hdist(single_cbs, j):
    single_cbs_list = list(single_cbs)
    if single_cbs_list[j-1] == '0':
        single_cbs_list[j-1] = '1'
        output_cbs = ''.join(single_cbs_list)
    else:
        single_cbs_list[j-1] = '0'
        output_cbs = ''.join(single_cbs_list)
    return output_cbs


def create_sigma_dict(sigma, Q_num):
    index_list = [i+1 for i in range(Q_num)]
    bit_list = [int(s) for s in list(sigma)]
    sigma_dict = dict(zip(index_list, bit_list))
    return sigma_dict


def create_tmp_sol_dict(GRB_dict):
    sol_dict = {}
    key_list = list(GRB_dict.keys())
    for key in key_list:
        sol_dict[key] = GRB_dict[key].x
    return sol_dict


def create_final_sol_dict(m, x, t, w, xi, Wbar, Tbar):
    var_class = ['Z', 'x', 't', 'w', 'xi', 'Wbar', 'Tbar']
    sol_dict = {}
    for V in var_class:
        if V == 'x': sol_dict[V] = create_tmp_sol_dict(x)
        if V == 't': sol_dict[V] = create_tmp_sol_dict(t)
        if V == 'w': sol_dict[V] = create_tmp_sol_dict(w)
        if V == 'xi': sol_dict[V] = create_tmp_sol_dict(xi)
        if V == 'Wbar': sol_dict[V] = create_tmp_sol_dict(Wbar)
        if V == 'Tbar': sol_dict[V] = create_tmp_sol_dict(Tbar)
        if V == 'Z': sol_dict[V] = m.ObjVal
    return sol_dict
