import numpy as np
import matplotlib.pyplot as plt

eps0=0.001
lambda_1=0.4
lambda_2=0.001
beta=5
k=0
x01=0
x02=1
x03=1
delta0=1
i=0

def func(x1,x2,x3):
    func = (x1**4)+17*(x2**4)+(x3**4)+24*(x1**2)*(x2**2)-8*x2*(x1**3)-32*x1*(x2**3)+6*(x2**2)*(x3**2)-4*x3*(x2**3)-4*x2*(x3**3)+3*(x1**2)+12*x1+2
    return func

def h1(x1_previous,x2_previous,x3_previous,eps01):
    x1_new=x1_previous+eps01
    func1=func(x1_new,x2_previous,x3_previous)
    func2=func(x1_previous,x2_previous,x3_previous)
    h1=(-1)*(func1-func2)/eps01
    return h1

def h2(x1_previous,x2_previous,x3_previous,eps01):
    x2_new=x2_previous+eps01
    h2=(-1)*(func(x1_previous,x2_new,x3_previous)-func(x1_previous,x2_previous,x3_previous))/eps01
    return h2

def h3(x1_previous,x2_previous,x3_previous,eps01):
    x3_new=x3_previous+eps01
    h3=(-1)*(func(x1_previous,x2_previous,x3_new)-func(x1_previous,x2_previous,x3_previous))/eps01
    return h3

def delta(x1_previous,x2_previous,x3_previous,h1,h2,h3,eps01):
    x1_new=x1_previous+eps01*beta*h1
    x2_new=x2_previous+eps01*beta*h2
    x3_new=x3_previous+eps01*beta*h3
    delta_new=func(x1_new,x2_new,x3_new)-func(x1_previous,x2_previous,x3_previous)
    return delta_new

'''
def loop01(x1,x2,x3,eps01,delta01):
    h1_new=h1(x1,x2,x3,eps01)
    h2_new=h2(x1,x2,x3,eps01)
    h3_new=h3(x1,x2,x3,eps01)
    delta_new=delta(x1,x2,x3,h1_new,h2_new,h3_new,eps01)
    return delta_new
'''

def psi_func(x1_previous,x2_previous,x3_previous,h1,h2,h3,delta,alpha,eps01):
    x1_new=x1_previous+alpha*h1
    x2_new=x2_previous+alpha*h2
    x3_new=x3_previous+alpha*h3
    psi=func(x1_new,x2_new,x3_new)-func(x1_previous,x2_previous,x3_previous)-alpha*delta*(1-lambda_1)/(beta*eps01)
    return psi

def phi_func(x1_previous,x2_previous,x3_previous,h1,h2,h3,delta,alpha,eps01):  
    x1_new=x1_previous+alpha*h1
    x2_new=x2_previous+alpha*h2
    x3_new=x3_previous+alpha*h3
    phi=func(x1_new,x2_new,x3_new)-func(x1_previous,x2_previous,x3_previous)-alpha*delta*lambda_1/(beta*eps01)
    return phi

def test_func(x1_previous,x2_previous,x3_previous,h1,h2,h3,p0,eps01):
    x1_new=x1_previous+p0*h1
    x2_new=x2_previous+p0*h2
    x3_new=x3_previous+p0*h3
    test=func(x1_new,x2_new,x3_new)-func(x1_previous,x2_previous,x3_previous)
    return test

eps=eps0
x1=x01
x2=x02
x3=x03
flag=0

while i<10:
    k=1000000000
    while True:
        p=1
        while delta0>=0.0:
            #delta0=loop01(x1,x2,x3,eps,delta0)
            h1_new=h1(x1,x2,x3,eps)
            h2_new=h2(x1,x2,x3,eps)
            h3_new=h3(x1,x2,x3,eps)
            delta0=delta(x1,x2,x3,h1_new,h2_new,h3_new,eps)
            if delta0>=0.0:
                eps=eps/2
        #print(delta0)
        #print(eps)

        mu=p

        while True:
            #number01=psi_func(x1,x2,x3,h1_new,h2_new,h3_new,delta0,mu,eps)
            while True:
                number01=psi_func(x1,x2,x3,h1_new,h2_new,h3_new,delta0,mu,eps)
                if number01>0.0:
                    break
                else:
                    mu=mu+1
            
            if number01==0.0:
                p=mu
                break

            number02=phi_func(x1,x2,x3,h1_new,h2_new,h3_new,delta0,mu,eps)

            if number02<=0.0:
                p=mu
                break

            a0=mu-p
            b0=mu

            while True:
                v0=(a0+b0)/2
                psi_v0=psi_func(x1,x2,x3,h1_new,h2_new,h3_new,delta0,v0,eps)
                phi_v0=phi_func(x1,x2,x3,h1_new,h2_new,h3_new,delta0,v0,eps)
                if psi_v0>=0.0 and phi_v0<=0.0:
                    p=v0
                    flag=1
                    break
                if psi_v0>0.0: 
                    b0=v0
                if psi_v0<=0.0:
                    a0=v0
            if flag==1:
                break
        
        k=test_func(x1,x2,x3,h1_new,h2_new,h3_new,p,eps)
        c=(p*lambda_1*delta0)/(beta*eps)
        if k<=(-1)*lambda_2*eps and k>=p*(1-lambda_1)*delta0 and k<=c:
            break
        else: eps=eps/2

        
    x1=x1+p*h1_new
    x2=x2+p*h2_new
    x3=x2+p*h3_new
    print("Ітерація "+str(i+1))
    print(x1,x2,x3)
    print(eps)


    i=i+1


c=func(x1,x2,x3)
print(c)
c=func(-2,-1,-1)
print(c)











'''
def test(delta1):
    i=0
    while delta1>1:
        delta1=delta1-3
        if delta1>1:
            i=i+1
    return i,delta1

i_new,delta_new=test(10)

print(i_new)
print(delta_new)
'''

'''
x1_new=x1+eps

h1=(-1*(func(x1_new,x2,x3)-func(x1,x2,x3)))/(eps**2)
'''
'''
h1=h1(x1,x2,x3,eps)
print(h1)
h2=h2(x1,x2,x3,eps)
print(h2)
h3=h3(x1,x2,x3,eps)
print(h3)
delta=delta(x1,x2,x3,h1,h2,h3,eps)
print(delta)
'''