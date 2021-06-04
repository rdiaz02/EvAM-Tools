#!/usr/bin/env python
from __future__ import division
import pdb
import matplotlib as mpl
import matplotlib
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import numpy as np
import pandas as pd
import argparse
import os
import string
import operator as op
from collections import OrderedDict
np.random.seed(0)

parser = argparse.ArgumentParser()
parser.add_argument("-f", required=True, default = None)
parser.add_argument("-transitions", required=False, default = None)
parser.add_argument("-end", required=False, default = "no", type=str)
parser.add_argument("-prob_type", required=False, default = "joint", type=str)
parser.add_argument("-data_type", required=False, default = "match-data", type=str)
parser.add_argument("-any_time", required=False, default = 0, type=int)
parser.add_argument("-weights", required=False, default="forwards_list-pord-edge-list-long-match-data.txt", type=str)
args = parser.parse_args()

args.transitions = "forwards_list-pord-trajectory-" + args.data_type + ".txt"
args.hypergraph_transitions = "forwards_list-pord-" + args.data_type + ".csv"

def LoadPi(in_file, N=100, rand=0):
    pis = []
    with open(in_file, 'r') as f:
        counter = 0
        for line_raw in f:
            if counter >= N:
                break
            line_raw2 = line_raw.rstrip("\r\n").split("\t")[1:]
            pi = []
            for line in line_raw2:
                pi.append([float(el) for el in line.split(" ")])
            pis.append(pi)
            counter += 1
    return pis
    
        
def MakeFeatureAdjacencyJoint(df, L, any_time=args.any_time, verbose=1):
    feature_transitions = [[0 for j in range(L+2)] for i in range(L+2)]
    count = 0
    for i in range(len(df)):
        traj = df.iloc[i,0]
        traj = [el for el in traj.split("-")[1:]]
        previous = 0
        if any_time == 0:
            for j, el in enumerate(traj):
                try:
                    el1 = int(el)+1
                except:
                    if str(el) != "end":
                        print el
                    el1 = L+1
                if previous == L+1:
                    print traj
                feature_transitions[previous][el1] += 1
                previous = el1
                count += 1
        else:
            temp = [0]
            for j, el in enumerate(traj):
                try:
                    el1 = int(el)+1
                except:
                    if str(el) != "end":
                        print el, type(el)
                    el1 = L+1
                for elp in temp:
                    feature_transitions[elp][el1] += 1
                temp.append(el1)
            count += 1
    for i, eli in enumerate(feature_transitions):
        for j, elj in enumerate(eli):
            feature_transitions[i][j] /= count
    return feature_transitions


def ProccessAdjacencyMatrix(f, transition_file, prob_type="joint", end = "yes", save="feature_transitions.csv"):
    pis = LoadPi(f, N=100)
    L = len(pis[0][0])

    df = pd.read_csv(transition_file, index_col=None, header = None)

    feature_transitions = MakeFeatureAdjacencyJoint(df, L)
    if end != "yes":
        feature_transitions = [el[:-1] for i, el in enumerate(feature_transitions[:-1])]
    if prob_type == "conditional":
        feature_transitions = [(el/np.sum(el) if np.sum(el) > 0 else el) for el in feature_transitions]
    feature_transitions = np.array(feature_transitions)
    
    np.savetxt(save, feature_transitions, delimiter=",")
    return feature_transitions

def GenotypeFromState(state):
    genes = np.array(list(string.ascii_uppercase))
    binary_state = reversed[list('{:0b}'.format(state))]
    genotype = genes[binary_state == 1]
    genotype = "".join(genotype) 
    return genotype

def CreateTransitionMatrix(file):
    data = pd.read_csv(file)
    data["from"] = data["from"].apply(GenotypeFromState)
    data["to"] = data["to"].apply(GenotypeFromState)


def main(args):
    ProccessAdjacencyMatrix(args.f, args.transitions, "conditional", end = args.end)
    transition_matrix = CreateTransitionMatrix(args.weight)

if __name__ == "__main__":
    main(args)

