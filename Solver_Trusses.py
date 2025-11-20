#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:34:19 2021

@author: kendrick
"""

import numpy as np

# compute unknown displacements 
def ComputeDisplacements(K, F, n_unknowns):
    # extract submatrix of unknowns
    K11 = K[0:n_unknowns,0:n_unknowns]
    F1 = F[0:n_unknowns]
    
    d = np.linalg.solve(K11,F1)
    
    return d

# postprocess the forces at known displacement nodes
def PostprocessReactions(K, d, F, n_unknowns, nodes):
    # These are computed net forces and do not
    # take into account external loads applied
    # at these nodes
    F = np.matmul(K[n_unknowns:,0:n_unknowns], d)
    
    # Postprocess the reactions
    for node in nodes:
        if node.xidx >= n_unknowns:
            node.AddReactionXForce(F[node.xidx-n_unknowns][0] - node.xforce_external)
        if node.yidx >= n_unknowns:
            node.AddReactionYForce(F[node.yidx-n_unknowns][0] - node.yforce_external)
        
    return F

# determine internal member loads

def ComputeMemberForces(bars):
    for bar in bars:
        E = bar.E
        A = bar.A
        L = bar.Length()
        lambdax, lambday = bar.LambdaTerms()
        
        u1, v1 = bar.init_node.xdisp, bar.init_node.ydisp
        u2, v2 = bar.end_node.xdisp, bar.end_node.ydisp
        
        delta = (-lambdax * u1 - lambday * v1 + lambdax * u2 + lambday * v2)
        F = ((E * A)/L) * delta
        
        bar.axial_load = F
        bar.is_computed = True
    pass


def ComputeNormalStresses(bars):
    
    for bar in bars:
        sigma = bar.axial_load/bar.A
        bar.normal_stress = sigma
        bar.is_computed = True
    pass

def ComputeBucklingLoad(bars):
    for bar in bars:
        L = bar.Length()
        Pcr = (np.pi**2 * bar.E * bar.Iu) / (144 * L**2)
        bar.buckling_load = Pcr
        bar.is_computed = True
    pass
