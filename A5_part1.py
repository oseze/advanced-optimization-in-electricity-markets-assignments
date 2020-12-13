# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 12:49:33 2019

@author: Pedro
"""

import gurobipy as gb

model_p1=gb.Model()
model_p2=gb.Model()
model_p3=gb.Model()
model_p4=gb.Model()
model_theta=gb.Model()


p1 = model_p1.addVar(lb=0, name='P1')

p2 = model_p2.addVar(lb=0, name='P2')

p3 = model_p3.addVar(lb=0, name='P3')

p4 = model_p4.addVar(lb=0, name='P4')

theta1 = model_theta.addVar(name='theta1')
theta2 = model_theta.addVar(name='theta2')

lam1 = 1
lam2 = 1
p_1 = 50
p_2 = 50
p_3 = 50
p_4 = 50
P1 = 0
P2 = 0
P = 5
lam1_old = 0.1
lam2_old = 0.1
old = 0
L1=200
L2=200
i = 1
gamma = 0.1
gamma1 = 0.2
a = 0
b = 1



model_p1.update()
model_p2.update()
model_p3.update()
model_p4.update()
model_theta.update()

#while (a+b*(2/gamma))>=0.000001:
while (abs((lam1-lam1_old))/abs(lam1_old))>= 0.000001 or (abs((lam2-lam2_old))/abs(lam2_old))>=0.000001:
    
    model_p1.setObjective(0.2*p1*p1+lam1*p1+(gamma/2)*(p1+(4*P-p_1)-400)*(p1+(4*P-p_1)-400),gb.GRB.MINIMIZE)
    c8=model_p1.addConstr(p1,gb.GRB.GREATER_EQUAL,0)
    
    model_p2.setObjective(0.1*p2*p2+lam1*p2+(gamma/2)*(p2+(4*P-p_2)-400)*(p2+(4*P-p_2)-400),gb.GRB.MINIMIZE)
    c9=model_p2.addConstr(p2,gb.GRB.GREATER_EQUAL,0)
        
    model_p3.setObjective(0.2*p3*p3+lam2*p3+(gamma/2)*(p3+(4*P-p_3)-400)*(p3+(4*P-p_3)-400),gb.GRB.MINIMIZE)
    c10=model_p3.addConstr(p3,gb.GRB.GREATER_EQUAL,0)
    
    model_p4.setObjective(0.5*p4*p4+lam2*p4+(gamma/2)*(p4+(4*P-p_4)-400)*(p4+(4*P-p_4)-400),gb.GRB.MINIMIZE)
    c11=model_p4.addConstr(p4,gb.GRB.GREATER_EQUAL,0)
    
    model_theta.setObjective(10*(theta1-theta2)*(lam2-lam1)+(gamma/2)*(p_1+p_2-10*(theta1-theta2)-L1)*(p_1+p_2-10*(theta1-theta2)-L1)+(gamma1/2)*(p_3+p_4-10*(theta2-theta1)-L2)*(p_3+p_4-10*(theta2-theta1)-L2),gb.GRB.MINIMIZE)
    c5=model_theta.addConstr(10*(theta2-theta1),gb.GRB.GREATER_EQUAL,-50)
    c6=model_theta.addConstr(10*(theta2-theta1),gb.GRB.LESS_EQUAL,50)
    c7=model_theta.addConstr(theta1,gb.GRB.EQUAL,0)
    
    model_p1.optimize()
    model_p2.optimize()
    model_p3.optimize()
    model_p4.optimize()
    model_theta.optimize()
    
    a = (gamma/2)*(p_1+p_2+10*(theta2.x-theta1.x)-L1)**2
    b = (gamma/2)*(p_3+p_4+10*(theta1.x-theta2.x)-L2)**2
    p_1 = p1.x
    p_2 = p2.x
    p_3 = p3.x
    p_4 = p4.x
    old = 10*(theta2.x-theta1.x)
    P = (p_1+p_2+p_3+p_4)/4.0
    P1 = (p_1+p_2)/2.0
    P2 = (p_3+p_4)/2.0
    lam1_old = lam1
    lam2_old = lam2
    #lam1 = lam1+(1/(1+0.1*i))*((L1-p1.x-p2.x-10*(theta2.x-theta1.x))/abs((L1-p1.x-p2.x-10*(theta2.x-theta1.x))))
    lam1 = lam1+gamma*(-L1+p1.x+p2.x-10*(theta1.x-theta2.x))
    #lam2 = lam2+(1/(1+0.1*i))*((L2-p3.x-p4.x-10*(theta1.x-theta2.x))/abs((L2-p3.x-p4.x-10*(theta1.x-theta2.x))))
    lam2 = lam2+gamma*(-L2+p3.x+p4.x-10*(theta2.x-theta1.x))
    i = i + 1
    
    model_p1.update()
    model_p2.update()
    model_p3.update()
    model_p4.update()
    model_theta.update()
    
    if i == 250:
        break
   
    
model_p1.setObjective(0.2*p1*p1+lam1*p1+(gamma/2)*(p1+(4*P-p_1)-400)*(p1+(4*P-p_1)-400),gb.GRB.MINIMIZE)
c8=model_p1.addConstr(p1,gb.GRB.GREATER_EQUAL,0)
    
model_p2.setObjective(0.1*p2*p2+lam1*p2+(gamma/2)*(p2+(4*P-p_2)-400)*(p2+(4*P-p_2)-400),gb.GRB.MINIMIZE)
c9=model_p2.addConstr(p2,gb.GRB.GREATER_EQUAL,0)
        
model_p3.setObjective(0.2*p3*p3+lam2*p3+(gamma/2)*(p3+(4*P-p_3)-400)*(p3+(4*P-p_3)-400),gb.GRB.MINIMIZE)
c10=model_p3.addConstr(p3,gb.GRB.GREATER_EQUAL,0)
    
model_p4.setObjective(0.5*p4*p4+lam2*p4+(gamma/2)*(p4+(4*P-p_4)-400)*(p4+(4*P-p_4)-400),gb.GRB.MINIMIZE)
c11=model_p4.addConstr(p4,gb.GRB.GREATER_EQUAL,0)
    
    
model_theta.setObjective(10*(theta1-theta2)*(lam2-lam1)+(gamma/2)*(p_1+p_2-10*(theta1-theta2)-L1)*(p_1+p_2-10*(theta1-theta2)-L1)+(gamma1/2)*(p_3+p_4-10*(theta2-theta1)-L2)*(p_3+p_4-10*(theta2-theta1)-L2),gb.GRB.MINIMIZE)
c12=model_theta.addConstr(10*(theta2-theta1),gb.GRB.GREATER_EQUAL,-50)
c13=model_theta.addConstr(10*(theta2-theta1),gb.GRB.LESS_EQUAL,50)
c14=model_theta.addConstr(theta1,gb.GRB.EQUAL,0)
    
model_p1.optimize()
model_p2.optimize()
model_p3.optimize()
model_p4.optimize()
model_theta.optimize()

print ('total cost', 0.2*p1.x**2+0.1*p2.x**2+0.2*p3.x**2+0.5*p4.x**2 )
print ('theta 2 ', theta2.x, 'theta1 ', theta1.x )