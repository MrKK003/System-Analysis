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
        self.xo = np.matrix([[xo], [xo], [xo]])
        self.xro = np.matrix([[1], [1], [1]])
        self.t_range = np.arange(0, 20 + self.to, self.to)
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
        plt.plot(x, self.y2)
        plt.plot(x, self.err)
        plt.legend(['x', 'x reconstructed', 'error'], loc=1)

        plt.show()

    # Algorithm core, loads three y arrays into class
    def generateY(self):
        phi = np.matrix(self.Phi())
        gamma = np.matrix(self.Gamma(self.Phi()))
        q = np.matrix([[2*self.to], [2*self.to], [self.to]])
        y = []
        y2 = []
        err = []
        num_of_iters = len(self.t_range)
        x_prev = self.xo
        xr_prev = self.xro
        u = 1
        y.append(x_prev[0])
        y2.append(xr_prev[0])
        err.append(np.linalg.norm(xr_prev - x_prev))
        for k in range(1, num_of_iters):
            if self.ko == 1:
                u = 0
            elif self.ko == 2:
                u = 3*np.sin(self.to*0.4*np.pi*k)
            elif self.ko == 3:
                if k % 2 == 0:
                    u = 2# Even 
                else:
                    u = -2 # Odd

            x = np.dot(phi, x_prev) + np.dot(gamma, u)
            xr = np.dot(phi, xr_prev) + np.dot(gamma, u) + np.dot(q,(x_prev[0] - xr_prev[0]))
            x_prev = x
            xr_prev = xr

            err.append(np.linalg.norm(xr - x))
            y.append(np.asscalar(x[0]))
            y2.append(np.asscalar(xr[0]))

        self.y1 = y
        self.y2 = y2
        self.err = err

    # calculates Phi to use in main formula
    def Phi(self):
        phi = np.matrix([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        q = self.q
        for i in range(q):
            phi += np.linalg.matrix_power(self.A*self.to, i+1)/np.math.factorial((i+1))
        return phi
        
    # calculates Gamma to use in main formula, needs Phi
    def Gamma(self, phi):
        return np.dot(np.dot((phi - np.identity(3)), np.linalg.inv(self.A)),self.B)

a = DynamicModel( 1, 3, 1, 10, 0.05, 1, 2)

a.runModel()