# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 12:49:33 2019

@author: Pedro
"""

import gurobipy as gb

model1=gb.Model()
model2=gb.Model()
model3=gb.Model()


#lam1 = model_1.addVar(lb=0,name='lam1_DA')
#lam2 = model_1.addVar(lb=0,name='lam2_DA')
p1_1 = model1.addVar(lb=0, ub = 50,  name='P1_1')
p2_1 = model1.addVar(lb=0, ub = 110, name='P2_1')
p3_1 = model1.addVar(lb=0, ub = 100, name='P3_1')
p1_2 = model2.addVar(lb=0, ub = 50,  name='P1_1')
p2_2 = model2.addVar(lb=0, ub = 110, name='P2_1')
p3_2 = model2.addVar(lb=0, ub = 100, name='P3_1')
p1_3 = model3.addVar(lb=0, ub = 50,  name='P1_1')
p2_3 = model3.addVar(lb=0, ub = 110, name='P2_1')
p3_3 = model3.addVar(lb=0, ub = 100, name='P3_1')
pw_1 = model1.addVar(lb=0, ub = 50, name='PW_1')
pw_2 = model2.addVar(lb=0, ub = 50, name='PW_2')
pw_3 = model3.addVar(lb=0, ub = 50, name='PW_3')
p3_rt1 = model1.addVar(lb=-45, ub = 45, name='P_RT1')
p3_rt2 = model2.addVar(lb=-45, ub = 45, name='P_RT2')
p3_rt3 = model3.addVar(lb=-45, ub = 45, name='P_RT3')
LS1 = model1.addVar(lb=0, ub=200, name='LS1')
SP1 = model1.addVar(lb=0, ub=50, name='SP1')
LS2 = model2.addVar(lb=0, ub=200, name='LS2')
SP2 = model2.addVar(lb=0, ub=22, name='SP2')
LS3 = model3.addVar(lb=0, ub=200, name='LS3')
SP3 = model3.addVar(lb=0, ub=10, name='SP3')


lam1_1 = 0.5
lam1_2 = 0.5
lam1_3 = 0.5
lam2_1 = 0.5
lam2_2 = 0.5
lam2_3 = 0.5
lam3_1 = 0.5
lam3_2 = 0.5
lam3_3 = 0.5
lamW_1 = 0.5
lamW_2 = 0.5
lamW_3 = 0.5
lam1_1old = 0.1
lam1_2old =0.1
lam1_3old =0.1 
lam2_1old =0.1 
lam2_2old =0.1
lam2_3old =0.1
lam3_1old =0.1
lam3_2old =0.1 
lam3_3old =0.1
lamW_1old =0.1
lamW_2old =0.1
lamW_3old =0.1
P1_DA = 0.5
P2_DA = 0.5
P3_DA = 0.5
PW_DA = 0.5
p_mean = 1
p_1 = 0.5
p_2 = 0.5
p_3 = 0.5
p_w = 0.5
i = 1
gamma1 = 0.1
gamma2 = 0.15




model1.update()
model2.update()
model3.update()

while (abs((lam1_1-lam1_1old))/abs(lam1_1old))>= 0.0001 or (abs((lam1_2-lam1_2old))/abs(lam1_2old))>= 0.0001 or (abs((lam1_3-lam1_3old))/abs(lam1_3old))>= 0.0001 or (abs((lam2_1-lam2_1old))/abs(lam2_1old))>= 0.0001 or (abs((lam2_2-lam2_2old))/abs(lam2_2old))>= 0.0001 or (abs((lam2_3-lam2_3old))/abs(lam2_3old))>= 0.0001 or (abs((lam3_1-lam3_1old))/abs(lam3_1old))>= 0.0001 or (abs((lam3_2-lam3_2old))/abs(lam3_2old))>= 0.0001 or (abs((lam3_3-lam3_3old))/abs(lam3_3old))>= 0.0001 or (abs((lamW_1-lamW_1old))/abs(lamW_1old))>= 0.0001 or (abs((lamW_2-lamW_2old))/abs(lamW_2old))>= 0.0001 or (abs((lamW_3-lamW_3old))/abs(lamW_3old))>= 0.0001:
    
    model1.setObjective(10*p1_1+25*p2_1+35*p3_1+0.2*(35*p3_rt1+200*LS1)+lam1_1*p1_1+lam2_1*p2_1+lam3_1*p3_1+lamW_1*pw_1+(gamma1/2)*((p1_1-P1_DA)*(p1_1-P1_DA)+(p2_1-P2_DA)*(p2_1-P2_DA)+(p3_1-P3_DA)*(p3_1-P3_DA))+(gamma2/2)*(pw_1-PW_DA)*(pw_1-PW_DA),gb.GRB.MINIMIZE)
    c1=model1.addConstr(p3_1+p3_rt1,gb.GRB.GREATER_EQUAL,0)
    c2=model1.addConstr(p3_1+p3_rt1,gb.GRB.LESS_EQUAL,100)
    c3=model1.addConstr(p1_1+p2_1+p3_1+pw_1,gb.GRB.EQUAL,200)
    c4=model1.addConstr(p3_rt1+(50-pw_1-SP1)+LS1,gb.GRB.EQUAL,0)
    
    model2.setObjective(10*p1_2+25*p2_2+35*p3_2+0.5*(35*p3_rt2+200*LS2)+lam1_2*p1_2+lam2_2*p2_2+lam3_2*p3_2+lamW_2*pw_2+(gamma1/2)*((p1_2-P1_DA)*(p1_2-P1_DA)+(p2_2-P2_DA)*(p2_2-P2_DA)+(p3_2-P3_DA)*(p3_2-P3_DA))+(gamma2/2)*(pw_2-PW_DA)*(pw_2-PW_DA),gb.GRB.MINIMIZE)
    c5=model2.addConstr(p3_2+p3_rt2,gb.GRB.GREATER_EQUAL,0)
    c6=model2.addConstr(p3_2+p3_rt2,gb.GRB.LESS_EQUAL,100)
    c7=model2.addConstr(p3_rt2+(22-pw_2-SP2)+LS2,gb.GRB.EQUAL,0)
    c8=model2.addConstr(p1_2+p2_2+p3_2+pw_2,gb.GRB.EQUAL,200) 
    
    model3.setObjective(10*p1_3+25*p2_3+35*p3_3+0.3*(35*p3_rt3+200*LS3)+lam1_3*p1_3+lam2_3*p2_3+lam3_3*p3_3+lamW_3*pw_3+(gamma1/2)*((p1_3-P1_DA)*(p1_3-P1_DA)+(p2_3-P2_DA)*(p2_3-P2_DA)+(p3_3-P3_DA)*(p3_3-P3_DA))+(gamma2/2)*(pw_3-PW_DA)*(pw_3-PW_DA),gb.GRB.MINIMIZE)
    c9=model3.addConstr(p3_3+p3_rt3,gb.GRB.GREATER_EQUAL,0)
    c10=model3.addConstr(p3_3+p3_rt3,gb.GRB.LESS_EQUAL,100)
    c11=model3.addConstr(p3_rt3+(10-pw_3-SP3)+LS3,gb.GRB.EQUAL,0)
    c12=model3.addConstr(p1_3+p2_3+p3_3+pw_3,gb.GRB.EQUAL,200) 
    
    model1.optimize()
    model2.optimize()
    model3.optimize()
    
    P1_DA = (0.2*p1_1.x+0.5*p1_2.x+0.3*p1_3.x)/3.0
    P2_DA = (0.2*p2_1.x+0.5*p2_2.x+0.3*p2_3.x)/3.0
    P3_DA = (0.2*p3_1.x+0.5*p3_2.x+0.3*p3_3.x)/3.0
    PW_DA = (0.2*pw_1.x+0.5*pw_2.x+0.3*pw_3.x)/3.0
    lam1_1old = lam1_1 
    lam1_2old =lam1_2 
    lam1_3old =lam1_3 
    lam2_1old =lam2_1 
    lam2_2old =lam2_2
    lam2_3old =lam2_3 
    lam3_1old =lam3_1 
    lam3_2old =lam3_2 
    lam3_3old =lam3_3 
    lamW_1old =lamW_1
    lamW_2old =lamW_2 
    lamW_3old =lamW_3
    
    lam1_1 = lam1_1+gamma1*(p1_1.x-P1_DA)
    lam1_2 = lam1_2+gamma1*(p1_2.x-P1_DA)
    lam1_3 = lam1_3+gamma1*(p1_3.x-P1_DA)
    lam2_1 = lam2_1+gamma1*(p2_1.x-P2_DA)
    lam2_2 = lam2_2+gamma1*(p2_2.x-P2_DA)
    lam2_3 = lam2_3+gamma1*(p2_3.x-P2_DA)
    lam3_1 = lam3_1+gamma1*(p3_1.x-P3_DA)
    lam3_2 = lam3_2+gamma1*(p3_2.x-P3_DA)
    lam3_3 = lam3_3+gamma1*(p3_3.x-P3_DA)
    lamW_1 = lamW_1+gamma2*(pw_1.x-PW_DA)
    lamW_2 = lamW_2+gamma2*(pw_2.x-PW_DA)
    lamW_3 = lamW_3+gamma2*(pw_3.x-PW_DA)
    
    i = i + 1
    
    model1.update()
    model2.update()
    model3.update()
    
    if i == 250:
        break
   
    
model1.setObjective(10*p1_1+25*p2_1+35*p3_1+0.2*(35*p3_rt1+200*LS1)+lam1_1*p1_1+lam2_1*p2_1+lam3_1*p3_1+lamW_1*pw_1+(gamma1/2)*((p1_1-P1_DA)*(p1_1-P1_DA)+(p2_1-P2_DA)*(p2_1-P2_DA)+(p3_1-P3_DA)*(p3_1-P3_DA))+(gamma2/2)*(pw_1-PW_DA)*(pw_1-PW_DA),gb.GRB.MINIMIZE)
c1=model1.addConstr(p3_1+p3_rt1,gb.GRB.GREATER_EQUAL,0)
c2=model1.addConstr(p3_1+p3_rt1,gb.GRB.LESS_EQUAL,100)
c3=model1.addConstr(p1_1+p2_1+p3_1+pw_1,gb.GRB.EQUAL,200)
c4=model1.addConstr(p3_rt1+(50-pw_1-SP1)+LS1,gb.GRB.EQUAL,0)
    
model2.setObjective(10*p1_2+25*p2_2+35*p3_2+0.5*(35*p3_rt2+200*LS2)+lam1_2*p1_2+lam2_2*p2_2+lam3_2*p3_2+lamW_2*pw_2+(gamma1/2)*((p1_2-P1_DA)*(p1_2-P1_DA)+(p2_2-P2_DA)*(p2_2-P2_DA)+(p3_2-P3_DA)*(p3_2-P3_DA))+(gamma2/2)*(pw_2-PW_DA)*(pw_2-PW_DA),gb.GRB.MINIMIZE)
c5=model2.addConstr(p3_2+p3_rt2,gb.GRB.GREATER_EQUAL,0)
c6=model2.addConstr(p3_2+p3_rt2,gb.GRB.LESS_EQUAL,100)
c7=model2.addConstr(p3_rt2+(22-pw_2-SP2)+LS2,gb.GRB.EQUAL,0)
c8=model2.addConstr(p1_2+p2_2+p3_2+pw_2,gb.GRB.EQUAL,200) 
    
model3.setObjective(10*p1_3+25*p2_3+35*p3_3+0.3*(35*p3_rt3+200*LS3)+lam1_3*p1_3+lam2_3*p2_3+lam3_3*p3_3+lamW_3*pw_3+(gamma1/2)*((p1_3-P1_DA)*(p1_3-P1_DA)+(p2_3-P2_DA)*(p2_3-P2_DA)+(p3_3-P3_DA)*(p3_3-P3_DA))+(gamma2/2)*(pw_3-PW_DA)*(pw_3-PW_DA),gb.GRB.MINIMIZE)
c9=model3.addConstr(p3_3+p3_rt3,gb.GRB.GREATER_EQUAL,0)
c10=model3.addConstr(p3_3+p3_rt3,gb.GRB.LESS_EQUAL,100)
c11=model3.addConstr(p3_rt3+(10-pw_3-SP3)+LS3,gb.GRB.EQUAL,0)
c12=model3.addConstr(p1_3+p2_3+p3_3+pw_3,gb.GRB.EQUAL,200) 

model1.optimize()
model2.optimize()
model3.optimize()

for v in model1.getVars():
    print(v.varName, v.x)

for v in model2.getVars():
    print(v.varName, v.x)
    
for v in model3.getVars():
    print(v.varName, v.x)

