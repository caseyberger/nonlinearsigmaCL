import numpy as np

class fields:
    def __init__(self,Nx,Ny,gL,n_eq,n_MC,integer = False):
        self.Nx = Nx
        self.Ny = Ny
        self.gL = gL
        self.n_eq = n_eq
        self.n_MC = n_MC
        self.representation = np.zeros([Nx,Ny,3])
        self.newlattice = np.zeros([Nx,Ny,3])
        self.integer = integer
        
    def random_spin(self):
        if self.integer:
            s = [0,0,0]
            i = np.random.randint(3)
            s[i] = 1
            return s
        else:
            s = np.random.rand(3)
            norm = np.linalg.norm(s)
            return s/norm
        
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
    
    def pick_next_triangle(self,x,y,vertex):
        x0 = vertex[0]
        y0 = vertex[1]
        if x0 == self.plus_one(x,self.Nx) and y0 == self.plus_one(y,self.Ny):
            return [self.minus_one(x0,self.Nx),y0]
        elif x0 == x and y0==self.plus_one(y,self.Ny):
            return [self.minus_one(x0,self.Nx),y0]
        elif x0==self.minus_one(x,self.Nx) and y0==self.plus_one(y,self.Ny):
            return [x0,self.minus_one(y0,self.Ny)]
        elif x0==self.minus_one(x,self.Nx) and y0 == y:
            return [x0,self.minus_one(y0,self.Ny)]
        elif x0==self.minus_one(x,self.Nx) and y0==self.minus_one(y,self.Ny):
            return [self.plus_one(x0,self.Nx),y0]
        elif x0==x and y0==self.minus_one(y,self.Ny):
            return [self.plus_one(x0,self.Nx),y0]
        elif x0==self.plus_one(x,self.Nx) and y0 ==self.minus_one(y,self.Ny):
            return [x0,self.plus_one(y0,self.Ny)]
        else:
            print("invalid point ({},{})".format(x0,y0))       
    
    def rho_squared(self,vertex_1, vertex_2, vertex_3):
        phi1 = self.phi(vertex_1[0],vertex_1[1])
        phi2 = self.phi(vertex_2[0],vertex_2[1])
        phi3 = self.phi(vertex_3[0],vertex_3[1])
        rhosq = 2.*(1.+ np.dot(phi1,phi2))*(1.+ np.dot(phi2,phi3))*(1.+ np.dot(phi3,phi1))
        return rhosq
    
    def q_sum(self,vertex_1, vertex_2, vertex_3):
        phi1 = self.phi(vertex_1[0],vertex_1[1])
        phi2 = self.phi(vertex_2[0],vertex_2[1])
        phi3 = self.phi(vertex_3[0],vertex_3[1])
        q = np.dot(phi1,phi2) + np.dot(phi2,phi3) + np.dot(phi3,phi1) + 1.j*np.dot(phi1,np.cross(phi2,phi3))
        return q
        
    def Q_L(self,newlattice = False):
        Q_L = 0
        triangle_dict = dict()
        triangle = 1
        vertex_1 = [x,y]
        vertex_2 = [self.plus_one(x,self.Nx),self.plus_one(y,self.Ny)]
        vertex_3 = [self.plus_one(x,self.Nx),y]
        triangle_dict[triangle] = [vertex_1, vertex_2, vertex_3]
        triangle+=1
        while triangle <= 8:
            q = self.q_sum(vertex_1, vertex_2, vertex_3)
            Q_L += -1.j*(np.log(1.+q)-0.5*np.log(rhosquared))/(2.*np.pi)
            vertex_3 = vertex_2
            vertex_2 = self.pick_next_triangle(x,y,vertex_2)
            triangle_dict[triangle] = [vertex_1, vertex_2, vertex_3]
            triangle+=1
        return Q_L
    
    def A_L(self,newlattice = False):
        A = 0.
        if newlattice:
            for x in range(self.Nx):
                for y in range(self.Nx):
                    xp1 = self.plus_one(x,self.Nx)
                    yp1 = self.plus_one(y,self.Ny)
                    A += np.dot(self.phi(x,y,newlattice = True),self.phi(xp1,y,newlattice = True))
                    A += np.dot(self.phi(x,y,newlattice = True),self.phi(x,yp1,newlattice = True))
        else:
            for x in range(self.Nx):
                for y in range(self.Nx):
                    xp1 = self.plus_one(x,self.Nx)
                    yp1 = self.plus_one(y,self.Ny)
                    A += np.dot(self.phi(x,y),self.phi(xp1,y)) 
                    A += np.dot(self.phi(x,y),self.phi(x,yp1))
        return A/self.gL
         
    def S_L(self,newlattice = False):
        S = self.A_L(newlattice = newlattice)-1.j*self.Q_L(newlattice = newlattice)
        return S
    
    def site(self,x,y):
        return self.representation[x][y]
    
    def Boltzmann_weight(self,x,y):
        S = self.A_L(x,y) - 1.j*self.Q_L(x,y)
        return np.exp(-1.*S)
    
    def change_spin(self,x,y):
        temp = self.representation
        temp[x][y] = self.random_spin()
        self.newlattice = temp
        oldS = self.S_L()
        newS = self.A_L(newlattice = True) - 1.j*self.Q_L(newlattice = True)
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
                for y in range(self.Nx):
                    self.change_spin(x,y)  
                    
    def Metropolis(self):
        mc_configs = []
        mc_S = []
        for n in range(self.n_MC):
            for x in range(self.Nx):
                for y in range(self.Nx):
                    self.change_spin(x,y) 
            mc_configs.append(self.representation)#.flatten())
            mc_S.append(self.S_L())
            if n%100 == 0:
                print("MC iteration {} complete".format(n))
        return mc_configs, mc_S
                    
    def check_Q(self):
        Q_err = []
        for x in range(self.Nx):
            for y in range(self.Ny):
                QL = self.Q_L(x,y)
                if np.absolute(QL)>0.5:
                    Q_err.append(QL)
        if len(Q_err) == 0:
            print("All values of Q are between -1/2 and 1/2")
        else:
            print("Some values of Q fall outside of -1/2 and 1/2")
        return Q_err