import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys 

class DynamicModel:
    def __init__(self, a1, a2, b, q, t, ko, xo):
        self.changeData(a1, a2, b, q, t, ko, xo)
    
    #1,1,1,3,0.002,0,1
    #1,1,1,3,0.002,2500,1
    #1,1,1,3,0.002,1800,1

    # just loads new data into the model
    def changeData(self, a1, a2, b, q, t, ko, xo):
        self.a1 = a1
        self.a2 = a2
        self.b = b
        self.q = q
        self.to = t
        self.ko = self.check_k(ko)
        self.xo = np.array([xo, xo, xo])
        self.t_range = np.arange(0, 10 + self.to, self.to)
        self.A = np.matrix([[0, 1, 0], [0, 0, 1], [-1, -self.a1, -self.a2]])
        self.B = np.array([[0, 0, self.b]])
        self.C = np.array([1, 0, 0])

    def check_k(self, ko):
        new_k = 0

        if ko==1:
          new_k=0
          return new_k

        if ko==2:
          new_k=2500
          return new_k
    
        if ko==3:
          new_k=1800
          return new_k

    # Runs algorithm, plots results
    def runModel(self):
        sns.set()
        self.generateY()
        y = self.y1
        x = self.t_range
        axes = plt.axes()
        axes.set_xlim([-0.2, 10.2])
        axes.set_xticks(np.arange(0, 10.5, 0.5))
        plt.xlabel('t - time')
        plt.ylabel('y(t) - output process')
        plt.scatter(x, y)
        plt.scatter(x, self.y2, c='yellow')
        plt.scatter(x, self.y3, c='red')
        plt.legend(['x1', 'x2', 'x3'], loc=4)
        plt.show()

    # Algorithm core, loads three y arrays into class
    def generateY(self):

        file = open('output.txt', 'a')
        
        phi = np.squeeze(np.array(self.Phi()))
        gamma = np.squeeze(self.Gamma(self.Phi()))
        y = []
        y2 = []
        y3 = []
        k = len(self.t_range)
        u = 1
        ucounter = 0
        x = np.array([])
        x_prev = np.squeeze(self.xo)
        y.append(x_prev[0])
        y2.append(x_prev[1])
        y3.append(x_prev[2])

        sys.stdout = file

        for _ in range(1, k):
            x = x_prev.dot(phi) + gamma*u
            print(x)
            y.append(x.dot(self.C))
            y2.append(x[1])
            y3.append(x[2])
            x_prev = x
            ucounter += 1
            if self.ko > 0:
                if ucounter >= self.ko:
                    u *= -1
                    ucounter = 0

        file.close()

        self.y1 = y
        self.y2 = y2
        self.y3 = y3

            
        

    # calculates Phi to use in main formula
    def Phi(self):
        phi = np.matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        q = self.q
        for i in range(q):
            phi += ((self.A*self.to)**(i+1))/np.math.factorial((i+1))
        return phi

    # calculates Gamma to use in main formula, needs Phi
    def Gamma(self, phi):
        I = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        a =  (phi - I)*self.A**(-1)
        a1 = np.squeeze(np.array(a))
        b = np.squeeze(np.array(self.B))
        return a1.dot(b)
        

    # for debuging and testing purposes
    def printData(self):
        print('A ', self.A)
        print('B ', self.B)
        print('C ', self.C)
        print('T ', self.t_range)
        phi = self.Phi()
        print('Phi ', phi)
        print('Gamma ', self.Gamma(phi))

'''try:
    a1 = int(input('Input a1: '))
    a2 = int(input('Input a2: '))
    b = int(input('Input b: '))
    q = int(input('Input q: '))
    t = float(input('Input t: '))
    k0 = int(input('Input k0: '))
    x0 = float(input('Input x0: '))
except ValueError:
    print('Error')'''

try:
    k0 = int(input('Enter number 1-3 to choose the right case: '))
    
except ValueError:
    print('Error')

a = DynamicModel( 1, 1, 1, 3, 0.002, k0, 1)
a.printData()
a.runModel()
