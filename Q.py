#!/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys, os
import argparse

def uqi(x,y):
    """A Universal Image Quality Index, named Q."""
    x = np.array([float(i) for i in x])
    y = np.array([float(i) for i in y])
    cov = np.cov(x,y)
    a_var = cov[0][0]
    b_var = cov[1][1]
    ab_cov = cov[0][1]
    uqi = (       (4 * ab_cov * x.mean() * y.mean() ) 
        / ((a_var + b_var) * (x.mean() ** 2 + y.mean() ** 2))
        )
    return uqi

def rmse(x,y):
    x = np.array([float(i) for i in x])
    y = np.array([float(i) for i in y])
    return np.sqrt(((x - y) ** 2).mean())

def nrmse(x,y):
    """normalized root-mean-square deviation or error"""
    x = np.array([float(i) for i in x])
    y = np.array([float(i) for i in y])
    nrmse = np.sqrt(((x - y) ** 2).mean())/y.mean()
    return nrmse

def vs_fig(x,y,filename=None,title=None,xlabel=None,ylabel=None,M=None,m=None):
#    font = {'family': 'serif',
#            'color':  'darkred',
#            'weight': 'normal',
#            'size': 16,
#            }
    x = np.array([float(i) for i in x])
    y = np.array([float(i) for i in y])
    fig, ax = plt.subplots()
    ax.scatter(x, y, s=2,color='red',zorder=10)
    lims = [np.min([ax.get_xlim(), ax.get_ylim()]),
            np.max([ax.get_xlim(), ax.get_ylim()]),]
    if M:lims[1]= M
    if m:lims[0]= m
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    xynrmse = nrmse(x,y)
    xyrmse = rmse(x,y)
    xyuqi = uqi(x,y)
    #text = "NRMSE:%f\nRMSE:%f\nUQI:%f\n"%(xynrmse,xyrmse,xyuqi)
    text = "RMSE:%f\nUQI:%f\n"%(xyrmse,xyuqi)
    ax.text(.9, .25, text,
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if filename:
        fig.savefig(filename,dpi=300)
        with open(filename+'.log',"w") as f:
            f.write("fig:%s\n"%(title))
            f.write("NRMSE:%f  normalized root-mean-square deviation or error.\n"%(xynrmse))
            f.write(" RMSE:%f  root-mean-square deviation or error.\n"%(xyrmse))
            f.write("  UQI:%f  A Universal Image Quality Index.\n"%(xyuqi))
    else:
        return fig, ax

def getdata(input_filename):
    MAX_COLS = 20
    data = [ [] for i in range(MAX_COLS)]
    with open(input_filename) as f:
        for i in f:
            if i[0] == '#':continue
            dat = i.split()
            for idx, num in enumerate(dat):
                data[idx].append(num)
    return data

def argparser():
    usage = ("To draw a fig vs 2 cols.")
    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument('-f', help="input file.")
    import tempfile
    tf = tempfile.NamedTemporaryFile(prefix="fig.", suffix='.png',dir='.')
    tf = tf.name
    parser.add_argument('-o', default=tf, help="output fig name.")
    parser.add_argument('-x', type=int, help="col number of x, start for 1.")
    parser.add_argument('-y', type=int, help="col number of y, start for 1.")
    parser.add_argument('-t', "--title", help="title of fig.")
    parser.add_argument('-xl', "--xlabel", help="xlabel of fig.")
    parser.add_argument('-yl', "--ylabel", help="ylabel of fig.")
    parser.add_argument("-M", type=float, help="max of x&y")
    parser.add_argument("-m", type=float, help="min of x&y")

    return parser

if __name__=='__main__':
    parser = argparser()
    args = parser.parse_args()
    if len(sys.argv) < 4:
        parser.print_help()
        sys.exit()
    data = getdata(args.f)
    x = data[args.x-1]
    y = data[args.y-1]
    vs_fig(x, y,
           filename=args.o, 
           title=args.title, 
           xlabel=args.xlabel, 
           ylabel=args.ylabel,
           M=args.M,
           m=args.m
           )

