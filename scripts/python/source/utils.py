""" useful python functions """

def parse_parameters(par_file):
    """ Reads par_file and constructs a dictionary with pairs: parameter name - parameter value.
    Format of the lines in param_file: parameter = value # optional comments """ 
    
    dict_param = dict([(l.split("=")[0].strip(), l.split("#")[0].split("=")[1].strip()) for l in par_file.readlines() if "=" in l and l.strip()[0]!="#"])
    return dict_param

def update_parameters(args):
    """  Update parameters defined from command line """ 
    pnames = args[0::2]
    tuple_args = zip(pnames, args[1::2])
    file_param = filter(lambda x: x[0]== '-p', tuple_args)
    # read parameter file
    fname = file_param[0][1]
    param = parse_parameters(open(fname))
    # read additional parameters
    other_param = filter(lambda x: x[0]!= '-p', tuple_args)
    if len(other_param) > 0:
        # check additional parameter names
        new_param = []
        for t_arg in other_param:
            if '--' in t_arg[0] and t_arg[0].split('--')[1] in param.keys():
                new_param.append((t_arg[0].split('--')[1], t_arg[1]))
            else:
                print "parameter %s ignored" %t_arg[0]
        # update parameters
        new_dict = dict([t for t in new_param if t[0] in param.keys()])
        param.update(new_dict)
    else:
        new_dict = param
    return param, fname, new_dict


def is_sorted(l):
    return all(l[i] <= l[i+1] for i in xrange(len(l)-1))

def FDR(p_values, ordered=False):
    """ Returns a list of FDR-corrected pvalues """ 
    if not ordered:
        ordered = is_sorted(p_values)
    else:
        joined = [ (v,i) for i,v in enumerate(p_values) ]
        joined.sort()
        p_values = [ p[0] for p in joined ]
        indices = [ p[1] for p in joined ]

    m = len(p_values)
    if m == 0 or not p_values:
        return []

    tmp_fdrs = [p*m/(i+1.0) for (i, p) in enumerate(p_values)]
    fdrs = []
    cmin = tmp_fdrs[-1]
    for f in reversed(tmp_fdrs):
        cmin = min(f, cmin)
        fdrs.append( cmin)
    fdrs.reverse()

    if not ordered:
        new = [ None ] * len(fdrs)
        for v,i in zip(fdrs, indices):
            new[i] = v
        fdrs = new

    return fdrs


def main():
    pass

if __name__ == '__main__':
    main()
