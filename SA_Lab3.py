import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

x_final = np.matrix([
    [2],
    [0], 
    [0]
     ], dtype=float)
class DynamicModel:
    def __init__(self, a1, a2, b, q, t, ko, xo):
        self.changeData(a1, a2, b, q, t, ko, xo)
    
    # just loads new data into the model
    def changeData(self, a1, a2, b, q, t, ko, xo):
        self.a1 = a1
        self.a2 = a2
        self.b = b
        self.q = q
        self.to = t
        self.ko = ko
        self.xo = np.matrix([[2], [0], [0]])
        self.t_range = np.arange(0, 1000 + self.to, self.to)
        self.t_range = self.t_range[:self.ko]
        self.A = np.matrix([[0, 1, 0], [0, 0, 1], [-1, -self.a1, -self.a2]])
        self.B = np.matrix([[0], [0], [self.b]])
        self.C = np.matrix([[1], [0], [0]])

    # Runs algorithm, plots results
    def runModel(self):
        sns.set()
        self.generateY()
        y = self.y1
        x = self.t_range
        plt.xlabel('t - time')
        plt.ylabel('y(t) - output process')
                
        plt.plot(x, y)
        plt.plot(x, self.y2, c='red')
        plt.plot(x, self.y3, c='green')
        plt.plot(x, self.uk, c='purple')
        plt.legend(['x1', 'x2', 'x3', 'u(k)'], loc=2)
        plt.show()

    # Algorithm core, loads three y arrays into class
    def generateY(self):
        phi = np.matrix(self.Phi())
        gamma = np.matrix(self.Gamma(self.Phi()))
        l0 = self.l0()
        y = []
        y2 = []
        y3 = []
        uk = []
        num_of_iters = len(self.t_range)
        x_prev = self.xo

        y.append(x_prev[0])
        y2.append(x_prev[1])
        y3.append(x_prev[2])

        for k in range(1, num_of_iters):
            u = self.U(k-1, l0)
            x = np.dot(phi, x_prev) + np.dot(gamma, u)
            x_prev = x

            uk.append(np.asscalar(u))
            y.append(np.asscalar(x[0]))
            y2.append(np.asscalar(x[1]))
            y3.append(np.asscalar(x[2]))

        uk.append(self.U(num_of_iters-1, l0))
        self.uk = uk
        self.y1 = y
        self.y2 = y2
        self.y3 = y3

    # calculates Phi to use in main formula
    def Phi(self):
        phi = np.matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        q = self.q
        for i in range(q):
            phi += np.linalg.matrix_power(self.A*self.to, i+1)/np.math.factorial((i+1))
        return phi

    # calculates Phi inverse to use in the formula for G(k)
    def Phi_inv(self):
        phi = np.matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        q = self.q
        m = -1
        for i in range(q):
           phi += m*np.linalg.matrix_power(self.A*self.to, i+1)/np.math.factorial((i+1))
           m *= -1
        return phi       

    # calculates Gamma to use in main formula, needs Phi
    def Gamma(self, phi):
        return np.dot(np.dot((phi - np.identity(3)), np.linalg.inv(self.A)),self.B)

    # calculates U(k) to use in main formula
    def U(self, k, l0):
        return np.dot(np.transpose(self.G(k)), l0)    
    
    # calculates l0 to use in the formula for U(k)
    def l0(self):    
        Phi = self.Phi()
        Pi_prev = np.identity(3)
        for _ in range(self.ko):
            Pi = np.dot(Phi, Pi_prev)
            Pi_prev = Pi
        
        Sigma = np.zeros((3, 3))
        for j in range(1, self.ko):
            Sigma = Sigma + np.dot(self.G(j), self.G(j).transpose())

        L = np.linalg.inv(np.dot(Pi, Sigma))
        return np.dot(L, x_final)      
            
    # calculater G(k) which is used in formulas for U(k) and l0
    def G(self, j):
        gamma = self.Gamma(self.Phi())
        phi_inv = self.Phi_inv()

        P_prev = np.identity(3)
        P = np.identity(3)
        for _ in range(0,j):
            P = np.dot(phi_inv, P_prev)
            P_prev = P
        return np.dot(P,gamma)


#a = DynamicModel( 1, 3, 1, 10, 0.03, 103, 2)

a = DynamicModel( 1, 3, 1, 10, 0.05, 110, 0)

a.runModel()