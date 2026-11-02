#!/usr/bin/env python3
"""Exact rational verification of M=21 coverage using independently computed gap tables.
The CSVs contain (M,p,s,j,gap_n,gap_d) with coefficient u=j/10.
No expected coverage counts or exception list are used in the test.
"""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
import csv
ROOT=Path(__file__).resolve().parent
M=21
C=((M-1)**2)//4
intervals=set()
for name in ('odd_moduli_exact_gaps.csv','even_moduli_exact_gaps.csv'):
    with (ROOT/name).open() as fh:
        for row in csv.reader(fh):
            if not row: continue
            m,p,s,j,gn,gd=map(int,row)
            if m!=M or gn>s*gd:continue
            assert 4*s<p and gd>0
            u=F(j,10); g=F(gn,gd)
            radius=(s-g)/C
            intervals.add((p,s,u-radius,u+radius))
assert (33,8,F(199,700),F(221,700)) in intervals

def floor(q):return q.numerator//q.denominator
def ceil(q):return -floor(-q)
def covers(a,b,z):
    p,s,c,d=z
    T=F(p,2)
    A=F(s); B=T-A
    return ceil((b*c-a*B)/T)<=floor((b*d-a*A)/T)

# The broad p=33 certificate: a/2 + (11/350)b >= 33/2.
# All pairs failing it are finite and enumerated from the exact inequality.
exceptions=[]; total=0; covered=0
for a in range(1,33):
    b=1
    while F(a,2)+F(11*b,350)<F(33,2):
        if gcd(a,b)==1:
            total+=1
            if any(covers(a,b,z) for z in intervals):covered+=1
            else:exceptions.append((a,b))
        b+=1
print('Exact certified coefficient intervals:',len(intervals))
print('Reduced fractions outside broad inequality:',total)
print('Certified by interval-intersection:',covered)
print('Not certified by these intervals:',len(exceptions))
print('List:', ', '.join(f'{a}/{b}' for a,b in exceptions))
if len(exceptions)==0:print('UNIFORM M=21 CERTIFIED')
else:print('Uniform M=21 NOT yet certified by these tests')
