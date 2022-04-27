import numpy as np
import matplotlib.pyplot as plt
import random
import math
import sys

class DynamicM:
    def __init__(self,a1,a2,b,q,t,k0,x0,k):
        self.changeData(a1,a2,b,q,t,k0,x0,k)

    def changeData(self,a1,a2,b,q,t,k0,x0,k):
        self.a1=a1
        self.a2=a2
        self.b=b
        self.q=q
        self.T0=t
        self.k0=k0
        self.k=k
        self.x0 = np.matrix([[0],[0],[0]])
        self.t_range = np.arange(0,300 + self.T0, self.T0)
        self.A = np.matrix([[0,1,0],[0,0,1],[-1,-self.a1,-self.a2]])
        self.B = np.matrix([[0],[0],[self.b]])
        self.C = np.matrix([1,0,0])

    def CalcY(self):
        u =1
        F = np.identity(3)
        for index in range(1, q+1):
            F +=np.linalg.matrix_power(self.A * T0,index)/float(math.factorial(index))
        G = np.dot(np.dot(F-np.identity(3),np.linalg.inv(self.A)),B)
        x_next = np.dot(F,self.x0) + (G*u)
        vector_arr = []
        X_arr1 = []
        Y_arr1 = []
        X_arr2 = []
        Y_arr2 = []
        X_arr3 = []
        Y_arr3 = []
        file = open('output.txt','a')
        print("success")
        sys.stdout = file
        plt.xlabel('t - time')
        plt.ylabel('y(t) - output process')
        for n in range(1,self.k):
            x_prev = x_next
            x_next = np.dot(F,x_prev)+(G*u)
            Y_arr1.append(float(np.dot(self.C,x_next)))
            vector_arr.append(x_next)
            X_arr1.append(n*self.T0)
        self.y1 = Y_arr1
        plt.plot(X_arr1, Y_arr1)
        x_next = np.dot(F,self.x0) + (G*u)
        for n in range(1,self.k):
            if (n>=int(self.k/2)):
                u=-1
            x_prev = x_next
            x_next = np.dot(F,x_prev)+(G*u)
            Y_arr2.append(float(np.dot(self.C,x_next)))
            vector_arr.append(x_next)
            X_arr2.append(n*self.T0)
        self.y2 = Y_arr2
        plt.plot(X_arr2, Y_arr2)
        x_next = np.dot(F,self.x0) + (G*u)
        for n in range(1,self.k):
            if (n<int(self.k*1/3)):
                u=1
            if (n>=int(self.k/3) and n<int(self.k*2/3)):
                u=-1
            if (n>=int(self.k*2/3)):
                u = 1
            x_prev = x_next
            x_next = np.dot(F,x_prev)+(G*u)
            Y_arr3.append(float(np.dot(self.C,x_next)))
            vector_arr.append(x_next)
            X_arr3.append(n*self.T0)
        self.y3 = Y_arr3
        plt.plot(X_arr3, Y_arr3)
        plt.show()
        print("First observation: \n X: ",X_arr1,"\n Y:",Y_arr1, "\n")
        print("Second observation: \n X: ",X_arr2,"\n Y:",Y_arr2, "\n")
        print("Third observation: \n X: ",X_arr3,"\n Y:",Y_arr3, "\n")
        file.close()


    def Run(self):
        self.CalcY()
        y = self.y1
        x=self.t_range
        axes = plt.axes()
        axes.set_xlim([-0.2, 200.2])
        axes.set_xticks(np.arange(0, 200.5, 10))
        plt.xlabel('t - time')
        plt.ylabel('y(t) - output process')
        plt.plot(x, self.y1, c ='blue')
        plt.plot(x, self.y2,c='red')
        plt.plot(x, self.y3,c='green')
        plt.legend(['x1', 'x2', 'x3'], loc=4)
        plt.show()
n=3
m=1
l=1
a1 = 1 #random.randrange(1,11,1)
a2 = 3 #random.randrange(1,11,1)
b = 1
A = np.matrix([[0,1,0],[0,0,1],[-1,-a1,-a2]])
B = np.matrix([[0],[0],[b]])
C = np.matrix([1,0,0])
Trand = random.randrange(1, 1000,1)
T0 = 0.02 #(1/float(Trand))
q = 10 #random.randrange(2,11,1)
model = DynamicM(a1,a2,b,q,T0,0,0,10000)
model.CalcY()
