import numpy as np

class fields:
    def __init__(self,Nx,Ny,theta,gL,n_eq,n_MC, freq = 100):
        self.Nx = Nx
        self.Ny = Ny
        self.theta = theta
        self.gL = gL
        self.n_eq = n_eq
        self.n_MC = n_MC
        self.freq = freq
        self.representation = np.zeros([Nx,Ny,3])
        self.newlattice = np.zeros([Nx,Ny,3])
        
    def random_spin(self):
        r = np.random.rand(2)
        inclination = np.arccos(1. - 2. * r[0])
        azimuth = 2. * np.pi * r[1]
        x = np.sin(inclination) * np.cos(azimuth)
        y = np.sin(inclination) * np.sin(azimuth)
        z = np.cos(inclination)
        return [x,y,z]
        
    def initialize_lattice(self,checkpoint = False, config = None):
        if checkpoint:
            self.representation = config
        else:
            #L = self.representation
            for x in range(self.Nx):
                for y in range(self.Ny):
                    #L[x][y] = self.random_spin()
                    self.representation[x][y] = self.random_spin()
            #self.representation = L
        
    def phi(self,x,y,newlattice = False):
        if newlattice:
            return self.newlattice[x][y]
        else:
            return self.representation[x][y]
    
    def plus_one(self,x,Nx):
        x+=1
        return x%Nx

    def minus_one(self,x, Nx):
        x-=1
        x+=Nx
        return x%Nx
    
    
    def Q_triangle(self,vertex_1, vertex_2, vertex_3, use_arccos = False, newlattice = False):
        phi1 = self.phi(vertex_1[0],vertex_1[1], newlattice = newlattice)
        phi2 = self.phi(vertex_2[0],vertex_2[1], newlattice = newlattice)
        phi3 = self.phi(vertex_3[0],vertex_3[1], newlattice = newlattice)
        rho2 = 2.*(1.+ np.dot(phi1,phi2))*(1.+ np.dot(phi2,phi3))*(1.+ np.dot(phi3,phi1))
        expQ_Re = (1. +  np.dot(phi1,phi2) + np.dot(phi2,phi3) + np.dot(phi3,phi1))/np.sqrt(rho2)
        expQ_Im = np.dot(phi1, np.cross(phi2, phi3))/np.sqrt(rho2)
        if use_arccos:
            return np.arccos(expQ_Re)/(2.*np.pi)
        else:
            return np.arcsin(expQ_Im)/(2.*np.pi)
        
    def Q_L(self,use_arccos = False, newlattice = False, check_Q = False):
        Q_L = 0
        Q_errs = 0
        #sum over all triangles on the lattice (2 triangles per site) 
        for x in range(self.Nx):
            for y in range(self.Ny):
                #triangle 1 -- lower, CCW
                v1 = [x,y]
                v2 = [self.plus_one(x,self.Nx),y]
                v3 = [self.plus_one(x,self.Nx),self.plus_one(y,self.Ny)]
                Q_tri = self.Q_triangle(v1,v2,v3,use_arccos = use_arccos, newlattice = newlattice)
                if check_Q:
                    if np.absolute(Q_tri)>0.5:
                        print("Some values of Q fall outside of -1/2 and 1/2: {}".format(Q_tri))
                        Q_errs += 1
                Q_L += Q_tri
                
                #triangle 2 -- upper, CCW
                v1 = [x,y]
                v2 = [self.plus_one(x,self.Nx),self.plus_one(y,self.Ny)]
                v3 = [x,self.plus_one(y,self.Ny)]
                Q_tri = self.Q_triangle(v1,v2,v3,use_arccos = use_arccos, newlattice = newlattice)
                if check_Q:
                    if np.absolute(Q_tri)>0.5:
                        print("Some values of Q fall outside of -1/2 and 1/2: {}".format(Q_tri))
                        Q_errs += 1
                Q_L += Q_tri
            if check_Q:
                if (Q_errs > 0):
                    print("{} total values of Q fall outside of -1/2 and 1/2".format(Q_errs))
                else:
                    print("No values of Q fall outside of -1/2 and 1/2!")
                if Q_L%1 == 0:
                    print("Q_L is an integer!")
                else:
                    print("Q_L is a non-integer value: {}".format(Q_L))
                return Q_L, Q_errs
            else:
                return Q_L
    
    def A_L(self,newlattice = False):
        A = 0.
        for x in range(self.Nx):
            for y in range(self.Ny):
                xp1 = self.plus_one(x,self.Nx)
                yp1 = self.plus_one(y,self.Ny)
                A += np.dot(self.phi(x,y,newlattice = newlattice),self.phi(xp1,y,newlattice = newlattice))
                A += np.dot(self.phi(x,y,newlattice = newlattice),self.phi(x,yp1,newlattice = newlattice))
        return A/self.gL
         
    def S_L(self, use_arccos = False, newlattice = False):
        S = self.A_L(newlattice = newlattice)
        S += -1.j*self.theta*self.Q_L(use_arccos = use_arccos,newlattice = newlattice)
        return S
    
    def site(self,x,y):
        return self.representation[x][y]
    
    def change_spin(self,x,y):
        temp = self.representation
        temp[x][y] = self.random_spin()
        self.newlattice = temp
        oldS = self.A_L(newlattice = False) - 1.j*self.theta*self.Q_L(newlattice = False)
        newS = self.A_L(newlattice = True) - 1.j*self.theta*self.Q_L(newlattice = True)
        dS = newS-oldS
        if dS < 0:
            self.representation[x][y] = self.newlattice[x][y]
        else:
            r = np.random.rand()
            if r <= np.exp(-1.*dS):
                self.representation[x][y] = self.newlattice[x][y]
            else:
                pass
       
    def equilibrate(self):
        for n in range(self.n_eq):
            for x in range(self.Nx):
                for y in range(self.Ny):
                    self.change_spin(x,y)  
                    
    def Metropolis(self):
        mc_configs = []
        mc_S = []
        mc_Q = []
        for n in range(self.n_MC):
            for x in range(self.Nx):
                for y in range(self.Ny):
                    self.change_spin(x,y)
            if n%self.freq == 0:
                print("MC iteration {} complete".format(n))
                mc_configs.append(self.representation)#.flatten())
                A_L = self.A_L()
                Q_L = self.Q_L()
                S_L = A_L - 1.j*self.theta*Q_L
                mc_S.append(S_L)
                mc_Q.append(Q_L)
        return mc_configs, mc_S, mc_Q
                    
    def check_Q(self):
        print("Computing Q_L with arccos:")
        QLcos, Qcos_err = self.Q_L(use_arccos = False, newlattice = False, check_Q = True)
        print("Computing Q_L with arcsin:")
        QLsin, Qsin_err = self.Q_L(use_arccos = True, newlattice = False, check_Q = True)
        if QLcos!= QLsin:
            print("Q takes on different values if solved with arcsin v arccos. Arccos = {}, arcsin = {}".format(QLcos, QLsin))
        else:
            print("Q is the same whether computed with arccos or arcsin: {}".format(QL_sin))
        return QLcos, Qcos_err, QLsin, Qsin_err