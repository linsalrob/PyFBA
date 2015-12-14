"""

Get a list of neighboring reactions and the number of genomes that have
that reaction in it. This will allow us to gapfill our model with 
some predicted reactions

"""


import os
import sys
from servers.SAP import SAPserver

server=SAPserver()


def close_genera(genome_name):
    res = server.all_genomes({'-complete' : True})
    close = {}
    for g in res:
        if genome_name in res[g]:
            close[g] = res[g]

    return close


def rast_closest_genomes(rast_file):
    close = {}
    with open(rast_file, 'r') as rin:
        for l in rin:
            p = l.strip().split("\t")
            close[p[0]] = p[2]
    return close


def closest_functions(close):
    """
    Measure the functions and score them by closeness (i.e. number
    of genomes they are in)
    """
    fn = {}

    cl = []
    if type(close) == type({}):
        cl = close.keys()
    else:
        cl = close

    fids = server.all_features({ '-ids' : cl})
    all_fids = set()
    for g in fids:
        all_fids.update(fids[g])

    all_fids = list(all_fids)
    i=0
    roles={}
    while (i<len(all_fids)):
        fid_list = all_fids[i:i+50000]
        roletmp = server.ids_to_data({ '-ids' : fid_list, '-data' : ['function'] })
        roles.update(roletmp)
        i += 50000


    for genome in fids:
        for f in fids[genome]:
            r = roles[f][0][0]
            if r not in fn:
                fn[r] = set()
            fn[r].add(genome)

    number = {r : 1.0 * len(fn[r])/len(close) for r in fn }
    return number



if __name__ == '__main__':

    try:
        genera = sys.argv[1]
    except:
        sys.exit(sys.argv[0] + " <genera name or closest.genomes file>")

    if os.path.exists(genera):
        # this is a closest genomes file
        cg = rast_closest_genomes(genera)
    else:
        # this is presumed to be a genera
        sys.stderr.write("Checking seed for genera like " + genera)
        cg = close_genera(genera)


    num = closest_functions(cg)
    for n in num:
        print(n + "\t" + str(num[n]))



