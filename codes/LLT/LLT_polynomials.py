#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computing LLT Polynomials Using Lattice Models

Created on Sat Jun 20 21:22:24 2020

@author: Calvin Yost-Wolff & Sylvester W. Zhang
"""

from sympy import Symbol
from sympy.utilities.iterables import ordered_partitions

def Partitions(n,outer,inner):
    '''
    bounded partitions
    '''
    for p in ordered_partitions(n):
        p.reverse()
        if len(p) >= len(inner) and len(p) <= len(outer):
            bounded = True
            for i in range(len(p)):
                if p[i] > outer[i]:
                    bounded = False
                if i < len(inner):
                    if p[i]<inner[i]:
                        bounded = False
            if bounded:
                yield p


def vertex_weight(n,a,b,C,D,i):
    Q = Symbol('q')
    boole = True
    if n in D and a == 1:
        return 0
    if (n-len(C)) + len(D) + a + (1-b) != n+1:
        return 0
    for j in D:
        if 0<j<n:
            if j+1 in C:
                boole = boole
            elif 1==1:
                return 0
    for j in C:
        if 1<j<n+1:
            if j-1 in D:
                boole = boole
            elif 1==1:
                return 0
    wt = 1
    if n in D:
        wt = Symbol('z_'+str(i))
    r = len(D)
    if a ==1:
        wt = wt*(Q**r)
    if a ==0 and b==1:
        if n in D:
            wt = wt*(Q**(r-1))
    if a==0 and b==0:
        if n in D:
            wt = wt*(Q**(r-1))
    return wt

#computes lambd/nu
def compute_ribbon_schur(lambd, nu, n, zi, zmax):
    #n is the ribbon size
    #zi is the lowest z_i term you want (probably z0 or z1)
    #zmax is the largest z_i you want (maximal entry)
    if zmax == zi:
        return row_traverse(lambd, nu, n, zi)
    part = []
    lam1 = []
    nu1 = []
    for i in range(len(nu)):
        part.append(0)
        lam1.append(0)
        nu1.append(0)
    sizelam = 0
    for i in range(len(lambd)):
        lam1[i] = lambd[i] + 1
        sizelam = sizelam + lam1[i]
    sizenu = 0
    for i in range(len(nu)):
        nu1[i] = nu[i] + 1
        sizenu = sizenu + nu1[i]
    part = []
    for i in range(len(nu)):
        part.append(0)
    s = 0
    for i in range(1 + sizelam - sizenu):
        for prho in Partitions(i+sizenu, outer=lam1, inner = nu1):
            for j in range(len(prho)):
                part[j] = prho[j]- 1
            recurse = row_traverse(lambd, part, n, zi)
            if recurse == 0:
                s = s
            elif 1==1:
                s = s + recurse*compute_ribbon_schur(part, nu, n, zi + 1, zmax)
    return s

def row_traverse(lambd, nu, n, zi):
    #input lambd and nu to be the same length
    #nu is the layer on top, lambd is the layer on bottom
    start = 0
    if len(lambd) != len(nu):
        return 0
    if lambd == nu:
        return 1
    lamrho = []
    nurho = []
    for i in range(len(lambd)):
        lamrho.append(0)
        nurho.append(0)
    for i in range(len(lambd)):
        lamrho[i] = lambd[i] + len(lambd) - i-1
        if lamrho[i]>start:
            start = lamrho[i]
    for i in range(len(nu)):
        nurho[i] = nu[i] + len(nu) - i-1
        if lamrho[i]>start:
            start = lamrho[i]
    currentD = []
    currentC = []
    counter = start
    wt = 1
    a = 0
    b = 0
    while counter > -1:
        currentC = []
        invert = 0
        for i in currentD:
            if i<n:
                currentC.append(i+1)
        if counter in lamrho:
            invert = invert + 1
            a = 1
        elif 1==1:
            a = 0
        if counter in nurho:
            b = 1
        elif 1==1:
            b = 0
            invert = invert + 1
        if n in currentD:
            invert = invert + 1
        if invert == 0 or invert == 3:
            return 0
        if invert == 2:
            currentC.append(1)
        wt = wt*vertex_weight(n,a,b,currentC,currentD,zi)
        if wt == 0:
            return 0
        counter = counter - 1
        currentD = currentC
    return wt

def LLT(Lambda,Nu,n,m):
    if len(Lambda)>len(Nu):
        while len(Lambda)>len(Nu):
            Nu+=[0]
    else:
        while len(Lambda) < len(Nu):
            Lambda+=[0]
    llt=compute_ribbon_schur(Lambda,Nu,n,1,m)
    if llt!=0:
        return llt.expand()
    else:
        return 0

def Schur(Lambda,Nu,m):
    return LLT(Lambda,Nu,1,m)
    return LLT(Lambda,Nu,1,m)
