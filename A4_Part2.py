#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 09:37:21 2018

@author: osezeiyore
"""

import gurobipy as gb
import math
import pandas as pd
#import numpy as np 

#%% Compact Bender's decomposition

UBB=0
UB = math.inf
UB2 = math.inf
LB = -math.inf
eps = 1e-5 
a_down = -1e18

sum_mast = pd.DataFrame(columns = ['g1','g2','g3','WT','theta1','theta2','a','obj'])
sum_sub1 = pd.DataFrame(columns = ['UB','G31','SP1','SH1','g1_k','g2_k','g3_k','WT_k','theta1_k','theta2_k','thetaRT1','thetaRT2','rho1','rho2','rho3','rhoWT','rhot1','rhot2'])
sum_sub2 = pd.DataFrame(columns = ['UB','G31','SP1','SH1','g1_k','g2_k','g3_k','WT_k','theta1_k','theta2_k','thetaRT1','thetaRT2','rho1','rho2','rho3','rhoWT','rhot1','rhot2'])
sum_sub3 = pd.DataFrame(columns = ['UB','G31','SP1','SH1','g1_k','g2_k','g3_k','WT_k','theta1_k','theta2_k','thetaRT1','thetaRT2','rho1','rho2','rho3','rhoWT','rhot1','rhot2'])

k = 0

# initialize master problem
mast = gb.Model()
mast.setParam( 'OutputFlag', False )

g1_status = mast.addVar(vtype=gb.GRB.BINARY) 
g2_status = mast.addVar(vtype=gb.GRB.BINARY) 
g3_status = mast.addVar(vtype=gb.GRB.BINARY) 

g1 = mast.addVar(lb=0, ub=50)
g2 = mast.addVar(lb=0, ub=110)
g3 = mast.addVar(lb=0, ub=100)
WT = mast.addVar(lb=0, ub =50)
theta1 = mast.addVar()
theta2 = mast.addVar()
a = mast.addVar(lb = -gb.GRB.INFINITY)

obj_m = (10*g1+25*g2+35*g3) + a 
mast.setObjective(obj_m,gb.GRB.MINIMIZE)
mast.addConstr(g1+g2+1000*(theta2-theta1)==200)  
mast.addConstr(g3+WT-1000*(theta2-theta1)==0)  
mast.addConstr(g1 >= g1_status*20)
mast.addConstr(g1 <= g1_status*50)
mast.addConstr(g2 >= g2_status*20)
mast.addConstr(g2 <= g2_status*110)
mast.addConstr(g3 >= g3_status*10)
mast.addConstr(g3 <= g3_status*100)
mast.addConstr(a >= a_down)

mast.update()
g1_f=0
g2_f=0
g3_f=0
WT_f=0
theta1_f=0
theta2_f=0

while abs(UB-LB)>eps:
    #---------------------------  SUB PROBLEM 1  ------------------------------
    sub = gb.Model()
    sub.setParam( 'OutputFlag', False )
    
    # variables (i)
    g1_k = sub.addVar(lb = 0, ub=50)
    g2_k = sub.addVar(lb = 0, ub=110)
    g3_k = sub.addVar(lb=0, ub=100)
    WT_k = sub.addVar(lb = 0, ub=50)
    theta1_k = sub.addVar()
    theta2_k = sub.addVar()
    
    # varibles 
    G31 = sub.addVar(lb=-45, ub=45) 
    SP1 = sub.addVar(lb=0, ub=50)
    SH1 = sub.addVar(lb=0, ub=200)
    thetaRT1 = sub.addVar(lb=-gb.GRB.INFINITY, ub=gb.GRB.INFINITY)
    thetaRT2 = sub.addVar(lb=-gb.GRB.INFINITY, ub=gb.GRB.INFINITY)
    
    sub.setObjective(0.2*(35*G31+200*SH1),gb.GRB.MINIMIZE)
    
    # constraints  
    if (g3_f>0):
        sub.addConstr(g3_k+G31 >=10)
    if (g3_f==0):
        sub.addConstr(g3_k+G31 ==0)
    sub.addConstr(0 <= g3_k+G31)    
    sub.addConstr(g3_k+G31 <= 100, name='pi')
    sub.addConstr(SP1 <= 50, name='mu')
    sub.addConstr(SH1 <= 200, name='o')
    sub.addConstr(g3_k+G31+WT_k<= 200)
    sub.addConstr(G31+(50-WT_k-SP1)-1000*((thetaRT2-thetaRT1)-(theta2_k-theta1_k))==0)
    sub.addConstr(SH1+1000*((thetaRT2-thetaRT1)-(theta2_k-theta1_k))==0)    
    sub.addConstr(g1_k == g1_f, name='rho1')
    sub.addConstr(g2_k == g2_f, name='rho2')
    sub.addConstr(g3_k == g3_f, name='rho3')
    sub.addConstr(WT_k == WT_f, name='rhoWT')
    sub.addConstr(theta1_k == theta1_f, name='rhot1')
    sub.addConstr(theta2_k == theta2_f, name='rhot2')
    sub.addConstr(1000*(thetaRT2-thetaRT1) <=200)
    sub.addConstr(1000*(thetaRT2-thetaRT1) >=-200)
    sub.update()
    sub.optimize()
    UB1 = sub.objVal
    sum_sub1.loc[k,'UB'] = UB1
    sum_sub1.loc[k,'G31'] = G31.x
    sum_sub1.loc[k,'SP1'] = SP1.x
    sum_sub1.loc[k,'SH1'] = SH1.x
    sum_sub1.loc[k,'g1_k'] = g1_k.x
    sum_sub1.loc[k,'g2_k'] = g2_k.x
    sum_sub1.loc[k,'g3_k'] = g3_k.x
    sum_sub1.loc[k,'WT_k'] = WT_k.x
    sum_sub1.loc[k,'theta1_k'] = theta1_k.x
    sum_sub1.loc[k,'theta2_k'] = theta2_k.x
    sum_sub1.loc[k,'thetaRT1'] = thetaRT1.x
    sum_sub1.loc[k,'thetaRT2'] = thetaRT2.x
    sum_sub1.loc[k,'rho1'] = sub.getConstrByName('rho1').Pi
    sum_sub1.loc[k,'rho2'] = sub.getConstrByName('rho2').Pi
    sum_sub1.loc[k,'rho3'] = sub.getConstrByName('rho3').Pi
    sum_sub1.loc[k,'rhoWT'] = sub.getConstrByName('rhoWT').Pi
    sum_sub1.loc[k,'rhot1'] = sub.getConstrByName('rhot1').Pi
    sum_sub1.loc[k,'rhot2'] = sub.getConstrByName('rhot2').Pi
    del sub
    
    #---------------------------  SUB PROBLEM 2  ------------------------------
  
    sub2 = gb.Model()
    sub2.setParam( 'OutputFlag', False )
    
     # variables (i)
    g1_k = sub2.addVar(lb = 0, ub=50)
    g2_k = sub2.addVar(lb = 0, ub=110)
    g3_k = sub2.addVar(lb=0, ub=100)
    WT_k = sub2.addVar(lb = 0, ub=50)
    theta1_k = sub2.addVar()
    theta2_k = sub2.addVar()
    
    # varibles 
    G31 =sub2.addVar(lb=-45, ub=45) 
    SP1 = sub2.addVar(lb=0, ub=50)
    SH1 = sub2.addVar(lb=0, ub=200)
    thetaRT1 = sub2.addVar()
    thetaRT2 = sub2.addVar()

    sub2.setObjective(0.5*(35*G31+200*SH1),gb.GRB.MINIMIZE)     
 
    # constraints 
    
    if (g3_f>0):
        sub2.addConstr(g3_k+G31 >=10)
    if (g3_f==0):
        sub2.addConstr(g3_k+G31 ==0)
    sub2.addConstr(0 <= g3_k+G31)    
    sub2.addConstr(g3_k+G31 <= 100, name='pi')
    sub2.addConstr(SP1 <= 22, name='mu')
    sub2.addConstr(SH1 <= 200, name='o')
    sub2.addConstr(g3_k+G31+WT_k<= 200)
    sub2.addConstr(G31+(22-WT_k-SP1)-1000*((thetaRT2-thetaRT1)-(theta2_k-theta1_k))==0)
    sub2.addConstr(SH1+1000*((thetaRT2-thetaRT1)-(theta2_k-theta1_k))==0)
    sub2.addConstr(g1_k == g1_f, name='rho1')
    sub2.addConstr(g2_k == g2_f, name='rho2')
    sub2.addConstr(g3_k == g3_f, name='rho3')
    sub2.addConstr(WT_k == WT_f, name='rhoWT')
    sub2.addConstr(theta1_k == theta1_f, name='rhot1')
    sub2.addConstr(theta2_k == theta2_f, name='rhot2')
    sub2.addConstr(1000*(thetaRT2-thetaRT1) <=200)
    sub2.addConstr(1000*(thetaRT2-thetaRT1) >=-200)
    sub2.update()
    sub2.optimize()
    UB2 = sub2.objVal
    sum_sub2.loc[k,'UB'] = UB2
    sum_sub2.loc[k,'G31'] = G31.x
    sum_sub2.loc[k,'SP1'] = SP1.x
    sum_sub2.loc[k,'SH1'] = SH1.x
    sum_sub2.loc[k,'g1_k'] = g1_k.x
    sum_sub2.loc[k,'g2_k'] = g2_k.x
    sum_sub2.loc[k,'g3_k'] = g3_k.x
    sum_sub2.loc[k,'WT_k'] = WT_k.x
    sum_sub2.loc[k,'theta1_k'] = theta1_k.x
    sum_sub2.loc[k,'theta2_k'] = theta2_k.x
    sum_sub2.loc[k,'thetaRT1'] = thetaRT1.x
    sum_sub2.loc[k,'thetaRT2'] = thetaRT2.x
    sum_sub2.loc[k,'rho1']  = sub2.getConstrByName('rho1').Pi
    sum_sub2.loc[k,'rho2']  = sub2.getConstrByName('rho2').Pi
    sum_sub2.loc[k,'rho3']  = sub2.getConstrByName('rho3').Pi
    sum_sub2.loc[k,'rhoWT'] = sub2.getConstrByName('rhoWT').Pi
    sum_sub2.loc[k,'rhot1'] = sub2.getConstrByName('rhot1').Pi
    sum_sub2.loc[k,'rhot2'] = sub2.getConstrByName('rhot2').Pi
    del sub2
      
    #---------------------------  SUB PROBLEM 3  ------------------------------
    
    sub3 = gb.Model()
    sub3.setParam( 'OutputFlag', False )
    
    # variables (i)
    g1_k = sub3.addVar(lb = 0, ub=50)
    g2_k = sub3.addVar(lb = 0, ub=110)
    g3_k = sub3.addVar(lb=0, ub=100)
    WT_k = sub3.addVar(lb = 0, ub=50)
    theta1_k = sub3.addVar()
    theta2_k = sub3.addVar()
    
    # varibles 
    G31 =sub3.addVar(lb=-45, ub=45) 
    SP1 = sub3.addVar(lb=0, ub=50)
    SH1 = sub3.addVar(lb=0, ub=200)
    thetaRT1 = sub3.addVar()
    thetaRT2 = sub3.addVar()

    sub3.setObjective(0.3*(35*G31+200*SH1),gb.GRB.MINIMIZE)    
 
    # constraints 
    if (g3_f>0):
        sub3.addConstr(g3_k+G31 >=10)
    if (g3_f==0):
        sub3.addConstr(g3_k+G31 ==0)
    sub3.addConstr( 0 <= g3_k+G31)    
    sub3.addConstr(g3_k+G31 <= 100, name='pi')
    sub3.addConstr(SP1 <= 10, name='mu')
    sub3.addConstr(SH1 <= 200, name='o')
    sub3.addConstr(G31+(10-WT_k-SP1)-1000*((thetaRT2-thetaRT1)-(theta2_k-theta1_k))==0)
    sub3.addConstr(SH1+1000*((thetaRT2-thetaRT1)-(theta2_k-theta1_k))==0)
    sub3.addConstr(g3_k+G31+WT_k<= 200)
    sub3.addConstr(g1_k == g1_f, name='rho1')
    sub3.addConstr(g2_k == g2_f, name='rho2')
    sub3.addConstr(g3_k == g3_f, name='rho3')
    sub3.addConstr(WT_k == WT_f, name='rhoWT')
    sub3.addConstr(theta1_k == theta1_f, name='rhot1')
    sub3.addConstr(theta2_k == theta2_f, name='rhot2')
    sub3.addConstr(1000*(thetaRT2-thetaRT1) <=200)
    sub3.addConstr(1000*(thetaRT2-thetaRT1) >=-200)
    sub3.update()
    sub3.optimize()
    UB3 = sub3.objVal
    sum_sub3.loc[k,'UB'] = UB3
    sum_sub3.loc[k,'G31'] = G31.x
    sum_sub3.loc[k,'SP1'] = SP1.x
    sum_sub3.loc[k,'SH1'] = SH1.x
    sum_sub3.loc[k,'g1_k'] = g1_k.x
    sum_sub3.loc[k,'g2_k'] = g2_k.x
    sum_sub3.loc[k,'g3_k'] = g3_k.x
    sum_sub3.loc[k,'WT_k'] = WT_k.x
    sum_sub3.loc[k,'theta1_k'] = theta1_k.x
    sum_sub3.loc[k,'theta2_k'] = theta2_k.x
    sum_sub3.loc[k,'thetaRT1'] = thetaRT1.x
    sum_sub3.loc[k,'thetaRT2'] = thetaRT2.x
    sum_sub3.loc[k,'rho1'] = sub3.getConstrByName('rho1').Pi
    sum_sub3.loc[k,'rho2'] = sub3.getConstrByName('rho2').Pi
    sum_sub3.loc[k,'rho3'] = sub3.getConstrByName('rho3').Pi
    sum_sub3.loc[k,'rhoWT'] = sub3.getConstrByName('rhoWT').Pi
    sum_sub3.loc[k,'rhot1'] = sub3.getConstrByName('rhot1').Pi
    sum_sub3.loc[k,'rhot2'] = sub3.getConstrByName('rhot2').Pi
    del sub3
        
    #--------------------------- MASTER PROBLEM  ------------------------------


    mast.addConstr(a >=0.2*(35*sum_sub1.loc[k,'G31'] + 200*sum_sub1.loc[k,'SH1'])+ 0.5*(35*sum_sub2.loc[k,'G31']\
                        + 200*sum_sub2.loc[k,'SH1'])+ 0.3*(35*sum_sub3.loc[k,'G31'] + 200*sum_sub3.loc[k,'SH1'])\
                        + (sum_sub1.loc[k,'rho1']+sum_sub2.loc[k,'rho1']+sum_sub3.loc[k,'rho1'])*(g1-(sum_sub1.loc[k,'g1_k']+sum_sub2.loc[k,'g1_k']+sum_sub3.loc[k,'g1_k'])*(1/3))\
                        + (sum_sub1.loc[k,'rho2']+sum_sub2.loc[k,'rho2']+sum_sub3.loc[k,'rho2'])*(g2-(sum_sub1.loc[k,'g2_k']+sum_sub2.loc[k,'g2_k']+sum_sub3.loc[k,'g2_k'])*(1/3))\
                        + (sum_sub1.loc[k,'rho3']+sum_sub2.loc[k,'rho3']+sum_sub3.loc[k,'rho3'])*(g3-(sum_sub1.loc[k,'g3_k']+sum_sub2.loc[k,'g3_k']+sum_sub3.loc[k,'g3_k'])*(1/3))\
                        + (sum_sub1.loc[k,'rhoWT']+sum_sub2.loc[k,'rhoWT']+sum_sub3.loc[k,'rhoWT'])*(WT-(sum_sub1.loc[k,'WT_k']+sum_sub2.loc[k,'WT_k']+sum_sub3.loc[k,'WT_k'])*(1/3))\
                        + (sum_sub1.loc[k,'rhot1']+sum_sub2.loc[k,'rhot1']+sum_sub3.loc[k,'rhot1'])*(theta1-(sum_sub1.loc[k,'theta1_k']+sum_sub2.loc[k,'theta1_k']+sum_sub3.loc[k,'theta1_k'])*(1/3))\
                        + (sum_sub1.loc[k,'rhot2']+sum_sub2.loc[k,'rhot2']+sum_sub3.loc[k,'rhot2'])*(theta2-(sum_sub1.loc[k,'theta2_k']+sum_sub2.loc[k,'theta2_k']+sum_sub3.loc[k,'theta2_k'])*(1/3)))



    mast.update()
    mast.optimize()

    g1_f = g1.x
    g2_f = g2.x
    g3_f = g3.x
    WT_f = WT.x
    theta1_f = theta1.x
    theta2_f = theta2.x
    LB = a.x
    UB=sum_sub1.loc[k,'UB']+sum_sub2.loc[k,'UB']+sum_sub3.loc[k,'UB']

    
    sum_mast.loc[k,'g1'] = g1_f
    sum_mast.loc[k,'g2'] = g2_f
    sum_mast.loc[k,'g3'] = g3_f
    sum_mast.loc[k,'WT'] = WT_f
    sum_mast.loc[k,'theta1'] = theta1_f
    sum_mast.loc[k,'theta2'] = theta2_f
    sum_mast.loc[k,'a'] = LB
    sum_mast.loc[k,'obj'] = mast.objVal
    sum_mast.loc[k,'UB']=UB


    k += 1


print(sum_sub1)   
print(sum_sub2)  
print(sum_sub3)   
print(sum_mast) 
