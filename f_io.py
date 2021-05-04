"""

########################################################################
# definition of functions to read, process & write files (txt,xyz,...) #
########################################################################

"""

import sys
import numpy as np
from scipy.constants import value

a0 = value('Bohr radius')*1e10 # (A)
Ha = value('Hartree energy in eV')

def Lorentz_2d(X,Y,x0,y0,Dx,Dy,w=1):  
    # X,Y = np.meshgrid(x,y)
    return w/(((X-x0)/Dx)**2+((Y-y0)/Dy)**2+1)

###############################################################################
#                               Read                                          # 
###############################################################################

def read(file,index=0,num=False):
    
    output = []
    with open(file, 'r') as infile:
        for line in infile:
             line = line.split()
             if line==[] or len(line)<(index+1):
                 continue
             else:
                 output.append(line[index])
    infile.close()    
    if num:
        tmp = [float(i) for i in output]
        output = np.asarray(tmp)
    return output

def read_xyz(file,idx=0,ext='xyz'):
    
    atm,X,Y,Z = [],[],[],[]
    with open(file+'.%s'%ext, 'r') as infile:
        for line in infile:
             line = line.split()
             if line==[] or line[0].islower():
                 continue
             elif len(line)>3:
                 atm.append(line[idx]); 
                 x = float(line[idx+1]); X.append(x)
                 y = float(line[idx+2]); Y.append(y)
                 z = float(line[idx+3]); Z.append(z)
    infile.close()
    X,Y,Z = np.asarray(X),np.asarray(Y),np.asarray(Z)
    return atm,X,Y,Z

def read_cif(file,unit='Ang'):
    
    i,idx = 0,0
    atm = []
    X0,Y0,Z0 = [],[],[]
    length,angles = [],[]
    with open(file+'.cif', 'r') as infile:
        for line in infile:
             i += 1
             line = line.split()
             if line==[]:
                 continue
             elif any([i for i in line[:1] if 'cell_length' in i]):
                 length.append(float(line[1]))
             elif any([i for i in line[:1] if 'cell_angle' in i]):
                 angles.append(float(line[1]))
             elif line==['loop_']:
                 idx = i;
             elif line==['_atom_site_fract_x']:
                 idx = i-idx-1
             elif any(i in ['_','#'] for i in line[0]): # or len(line)<5 
                 continue
             elif len(line)>4:
                 atmi = ''.join([i for i in line[0] if not i.isdigit()])
                 atm.append(atmi); 
                 x = float(line[idx+0]); X0.append(x)
                 y = float(line[idx+1]); Y0.append(y)
                 z = float(line[idx+2]); Z0.append(z)
    infile.close()
    X0,Y0,Z0 = np.asarray(X0),np.asarray(Y0),np.asarray(Z0)
    
    a,b,c = length
    angles = np.asarray(angles)*np.pi/180
    alpha,beta,gamma = angles # rad

    dcos = np.cos(alpha)-np.cos(beta)*np.cos(gamma)
    fcos = 2*np.prod(np.cos(angles))-np.sum(np.cos(angles)**2)
    Omega = c*np.sqrt(1+fcos)
    
    cell = np.zeros((3,3))
    cell[0,:] = np.array([a,0,0])
    cell[1,:] = np.array([b*np.cos(gamma),b*np.sin(gamma),0])
    cell[2,:] = np.array([c*np.cos(beta),c*dcos/np.sin(gamma),Omega/np.sin(gamma)])
    
    if unit=='Ang':
        X = np.dot(cell[:,0],np.array([X0,Y0,Z0]))
        Y = np.dot(cell[:,1],np.array([X0,Y0,Z0]))
        Z = np.dot(cell[:,2],np.array([X0,Y0,Z0]))
        return atm,X,Y,Z,cell
    elif unit=='frac':
        return atm,X0,Y0,Z0,cell
    else:
        sys.exit('   unit is Ang or frac')

def read_cube(file,Na,f_type='gaussian'):
    
    # gaussian : http://paulbourke.net/dataformats/cube/
    
    n,mesh = np.zeros(3,dtype=int),np.zeros((3,3))
    output = []; 
    with open(file+'.cube', 'r') as infile:
        if f_type=='gaussian':
            idx = 1;
            for line in infile:
                line = line.split(); 
                if line==[]:
                    continue
                elif (idx-6)>Na :
                    for line_i in line:
                        output.append(float(line_i))
                elif idx==3:
                    line_f = [float(i) for i in line[1:]]
                    x0 = np.asarray(line_f)*a0
                elif idx in [4,5,6]:
                    line_f = [float(i) for i in line[1:]]
                    n[idx-4],mesh[idx-4,:] = int(line[0]),line_f
                idx += 1
        else:
            sys.exit('   cube file not implemented')
    infile.close()
    output_arr = np.zeros(n)
    for i in range(n[0]):
        for j in range(n[1]):
            for k in range(n[2]):
                output_arr[i,j,k] = output[(i*n[1]+j)*n[2]+k]
    return output_arr,n,x0 

###############################################################################
#                               Print                                         # 
###############################################################################

def print_xyz(atm,x,y,z):
    for atmi,xi,yi,zi in zip(atm,x,y,z):
        str_xi = '{0:.10f}'.format(xi)
        str_yi = '{0:.10f}'.format(yi)
        str_zi = '{0:.10f}'.format(zi)
        print('   ',atmi,'   ',str_xi,'   ',str_yi,'   ',str_zi);

def print_md2(atm,x,y,z,q0=[],i0=1):
    print()
    if q0==[]:
        q0 = OpenMX.q0(None,atm=atm)
    for i,atmi in enumerate(atm):
        str_xi = '{0:.10f}'.format(x[i])
        str_yi = '{0:.10f}'.format(y[i])
        str_zi = '{0:.10f}'.format(z[i])
        str_atm = '   '+str(i+i0)+'   '+str(atmi);
        str_xyz = '   '+str_xi+'   '+str_yi+'   '+str_zi;
        str_q0  = '   '+str(q0[i])+' '+str(q0[i]);
        print(str_atm+str_xyz+str_q0);

###############################################################################
#                               Write                                         # 
###############################################################################

def write_xyz(file,atm,x,y,z):
    file.write('\n \n')
    for atmi,xi,yi,zi in zip(atm,x,y,z):
            file.write('   '+str(atmi)+'  '+str(xi)+'  '+str(yi)+'  '+str(zi)+'\n');    

def write_md2(file,atm,x,y,z,q0=[],i0=1):
    if q0==[]:
        q0 = OpenMX.q0(None,atm=atm)
    for i,atmi in enumerate(atm):
            str_xi = '{0:.4f}'.format(x[i])
            str_yi = '{0:.10f}'.format(y[i])
            str_zi = '{0:.10f}'.format(z[i])
            str_atm = '   '+str(i+i0)+'   '+str(atmi);
            str_xyz = '   '+str_xi+'   '+str_yi+'   '+str_zi;
            str_q0  = '   '+str(q0[i])+' '+str(q0[i]);
            file.write(str_atm+str_xyz+str_q0+'\n');
            
def write_cube(file,X,U,dim=3):
    if dim==3:
        x,y,z = X
        for xi,yi,zi,ui in zip(x,y,z,U):
            file.write('  '+str(xi)+'  '+str(yi)+'  '+str(zi)+'  '+str(ui)+'\n')
    elif dim==2:
        x,y = X
        for xi,yi,ui in zip(x,y,U):
            file.write('  '+str(xi)+'  '+str(yi)+'  '+str(ui)+'\n')        

###############################################################################
#                               OpenMX                                        # 
###############################################################################

class OpenMX:
    
    def __init__(self,name,path='.',ext='out',read=True):
        
        # file
        self.name = name
        self.path = '%s/'%path
        self.ext = ext
        # system
        self.Ns = 0 # Species number
        self.species = []
        self.Na = 0 # Atoms number
        self.atm = []
        self.bs = [] # basis set
        self.acell = np.zeros((3,3))
        self.ngrid = np.zeros(3,dtype=int) # number of grid
        self.gcell = np.zeros((3,3))
        # self-consistent cycle
        self.solver = ''
        self.Nk = np.zeros(3,dtype=int)
        # NEGF 
        self.Nk_gf = np.zeros(2,dtype=int)
        self.NEt = 0
        self.Et = np.zeros(0)
        self.Nkt = np.zeros(2,dtype=int)
        # energy
        self.Utot = 0
        self.Ekin = 0
        self.mu = 0
        self.mu_ll = 0 # mu left lead
        self.mu_rl = 0 # mu right lead
        # orbitals
        self.s = []
        self.p = []
        self.d = []
        # orbitals unfolding
        self.orb_id = []
        self.Ebounds = np.zeros(2)
        self.Nk_orb = 0
        self.k_orb = np.zeros(0)
        self.Nkpoints = 0 # Nk in path definition
        self.Nkpath_orb = np.zeros(0,dtype=int) # Nk for each path
        self.ikpath_orb = np.zeros(0,dtype=int) # index of points in path definition 
        # Mulliken population
        self.MC = np.zeros(0) # Na
        # geometry
        self.coord = np.zeros((0,3)) # Na,3
        self.force = np.zeros((0,3)) # Na,3
        # band structure
        self.nbds = 0
        self.Nkpath = 0
        self.kpath = np.zeros(0,dtype=int) # 2,Nkapth
        
        # read
        if read:
            self.read()
            if (self.acell==0).all():
                self.acell = self.gcell*self.ngrid*a0
    
    def read(self):
        
        read_Ek = False
        read_orb = False
        infile = open(self.path+self.name+'.%s'%self.ext,'r')    
        for line in infile:
            line = line.split()
            
            # self-consistent cycle
            if line[:1]==['scf.EigenvalueSolver']:
                self.solver = line[1]
            elif line[:1]==['scf.Kgrid']:
                line = [int(l) for l in line[1:4]]
                self.Nk = np.asarray(line)
            elif line[:1]==['NEGF.scf.Kgrid']:
                line = [int(l) for l in line[1:3]]
                self.Nk_negf = np.asarray(line)
                
            # transmission
            elif line[:1]==['NEGF.tran.energydiv']:
                self.NEt = int(line[1])
            elif line[:1]==['NEGF.tran.energyrange']:
                self.Et = np.linspace(float(line[1]),float(line[2]),self.NEt)
                self.nu = float(line[3])
            elif line[:1]==['NEGF.tran.Kgrid']:
                line = [int(l) for l in line[1:3]]
                self.Nkt = np.asarray(line)
                
            # system
            elif line[:1]==['Species.Number']:
                self.Ns = int(line[1])
            elif line[:1]==['<Definition.of.Atomic.Species']:
                for _,line in zip(range(self.Ns),infile):
                    line = line.split()
                    self.species.append(line[0])
                    id_bs = line[1].index('-')+1
                    self.bs.append(line[1][id_bs:id_bs+6])
            elif line[:1]==['Atoms.Number']:
                self.Na += int(line[1])
            elif line[:1]==['LeftLeadAtoms.Number']:
                self.Na += int(line[1])
            elif line[:1]==['RightLeadAtoms.Number']:
                self.Na += int(line[1])                
            elif line[:1]==['<Atoms.UnitVectors']:
                for i,line in zip(range(3),infile):
                    line = line.split()
                    line = [float(l) for l in line]
                    self.acell[i,:] = line; 
            elif line[:3]==['Num.','of','grids']:
                line = [int(l.replace(',','')) for l in line[9:]]
                self.ngrid = np.asarray(line)
            elif line[:3]==['Cell','vectors','(bohr)']:
                for i,line in zip(range(3),infile):
                    line = line.split()
                    line = [float(l.replace(',','')) for l in line[2:]]
                    self.gcell[i,:] = line
                    
            # energy
            elif line[:1]==['Utot.']:
                self.Utot = float(line[1])*Ha
            elif line[:1]==['Ukin.']:
                self.Ekin = float(line[1])*Ha                
            elif line[:2]==['Chemical','Potential']:
                self.mu = float(line[4])*Ha
            elif line[:5]==['Chemical','potential','of','left','lead']:
                self.mu_ll = float(line[6])*Ha
            elif line[:5]==['Chemical','potential','of','right','lead']:
                self.mu_rl = float(line[6])*Ha            
                
            # band structure
            elif line[::2]==['k1=','k2=','k3=']:
                read_Ek = not read_Ek
                next(infile)
            elif read_Ek:
                if line==[]:
                    read_Ek = not read_Ek
                    if self.solver=='NEGF':
                        nlines = ((np.prod(self.Nk_negf)+1)//2-1)*(4+self.nbds)
                    else:
                        nlines = ((np.prod(self.Nk)+1)//2-1)*(4+self.nbds)
                    for _ in range(nlines):
                        next(infile)
                    continue
                self.nbds = int(line[0])
                
            # bands unfolding
            elif line[:1]==['Unfolding.LowerBound']:
                self.Ebounds[0] = float(line[1])
            elif line[:1]==['Unfolding.UpperBound']:
                self.Ebounds[1] = float(line[1])
            elif line[:1]==['Unfolding.Nkpoint']:
                self.Nkpoints = int(line[1])-1
            elif line==['Unfolding','calculation','for','band','structure']:
                read_orb = not read_orb
                orb_i = []
                for _ in range(11):
                        next(infile)
            elif read_orb:
                if line==[]:
                    read_orb = not read_orb
                    self.orb_id.append(orb_i)
                    self.orb_id = self.orb_id[1:]
                    continue
                elif len(line)==5:
                    self.orb_id.append(orb_i);
                    orb_i = []
                orb_i.append(int(line[0]));
            elif line[:3]==['For','each','path:']:
                line = [int(l) for l in line[4:4+self.Nkpoints+1]]
                self.Nkpath_orb = np.asarray(line,dtype=int)
                self.Nk_orb = np.sum(self.Nkpath_orb)
                line = [np.sum(line[:i]) for i,_ in enumerate(line)]
                self.ikpath_orb = np.asarray(line,dtype=int)
            elif line==['ka','kb','kc']:
                self.k_orb = np.zeros((self.Nk_orb,3))
                for i,line in zip(range(self.Nk_orb),infile):
                    line = line.split()
                    line = [float(l) for l in line[-3:]]
                    self.k_orb[i,:] = line
                
            # Mulliken population
            elif line==['Up','spin','Down','spin','Sum','Diff']:
                self.MC = np.zeros((self.Na,3))
                for i,line in zip(range(self.Na),infile):
                    line = line.split()
                    self.atm.append(line[1])
                    line = [float(l) for l in line[2:5]]
                    self.MC[i,:] = line
                    
            # orbitals (from Mulliken)
            elif line==['Decomposed','Mulliken','populations']:
                j,s,p,d = 0,[],[],[]
                for atm_i in self.atm:
                    idx = self.species.index(atm_i)
                    tmp = self.bs[idx][1::2]
                    nbs = int(tmp)*10**(3-len(tmp))
                    a,b,c = str(nbs)
                    Ns,Np,Nd = int(a),int(b),int(c)
                    nlines = 2*Ns+4*Np+6*Nd+len(tmp)+3;
                    for i,line in zip(range(nlines),infile):
                        line = line.split()
                        if line[:1]==['s']:
                            s.append(j); j += 1
                        elif line[:1] in [['px'],['py'],['pz']]:
                            p.append(j); j += 1
                        elif line[:1] in [['d3z^2-r^2'],['dx^2-y^2'],['dxy'],['dxz'],['dyz']]:
                            d.append(j); j+= 1
                        if i==nlines-1:
                            self.s.append(s)
                            self.p.append(p)
                            self.d.append(d)
                            s,p,d = [],[],[]
                self.s = np.asarray(self.s)
                self.p = np.asarray(self.p)
                self.d = np.asarray(self.d)
                    
            # coordinates & forces
            elif line[:1]==['<coordinates.forces']:
                self.coord = np.zeros((self.Na,3))
                self.force = np.zeros((self.Na,3))
                next(infile)
                for i,line in zip(range(self.Na),infile):
                    line = line.split()
                    line = [float(l) for l in line[2:]]
                    self.coord[i,:],self.force[i,:] = line[:3],line[3:]
                    
        infile.close()
    
    def q0(self,atm=''):   
        
        if atm=='':
            atm = self.atm   
            
        q = []
        for atmi in atm:
            if atmi[:1]=='E':
                q.append(0.0) # Empty atoms
            elif atmi in ['H']:
                q.append(0.5)
            elif atmi in ['Al']:
                q.append(1.5)
            elif atmi in ['C']:
                q.append(2.0)
            elif atmi in ['N']:
                q.append(2.5)
            elif atmi in ['O','S','Se']:
                q.append(3.0)
            elif atmi in ['Cu']:
                q.append(5.5)
            elif atmi in ['Ti','W']:
                q.append(6.0)
            elif atmi in ['Mo','Cr']:
                q.append(7.0)
            elif atmi in ['Pd','Pt']:
                q.append(8.0)
            elif atmi in ['Au','Ag']:
                q.append(8.5)
        return np.asarray(q)   
        
    def Ek(self,ext='out'):
        
        # read_Ek
        ik = 0
        read_Ek = False
        nk = (np.prod(self.Nk)+1)//2
        k0 = np.zeros((3,nk))
        Ek0 = np.zeros((nk,self.nbds))
        infile = open(self.path+self.name+'.%s'%ext,'r')    
        for line in infile:
            line = line.split()
            if line[::2]==['k1=','k2=','k3=']:
                read_Ek = not read_Ek
                line = [float(l) for l in line[1::2]]
                k0[:,ik] = np.asarray(line)
                continue
            elif read_Ek:
                if line==[]:
                    continue
                elif line==['nbands']:
                    read_Ek = not read_Ek
                    ik += 1     
                    continue
                iE = int(line[0])
                Ek0[ik,iE-1] = float(line[1])*Ha 
                if iE==self.nbds:
                    read_Ek = not read_Ek
                    ik += 1
        
        kx,ikx = np.unique(k0[0],return_inverse=True)
        ky,iky = np.unique(k0[1],return_inverse=True)
        kz,ikz = np.unique(k0[2],return_inverse=True)
        nkx,nky,nkz = len(kx),len(ky),len(kz)
        
        Ek = np.ones((nkx,nky,nkz,self.nbds))*100 # eV
        for i,Ei in enumerate(Ek0):
            ix,iy,iz = ikx[i],iky[i],ikz[i]
            Ek[ix,iy,iz] = Ek0[i]
        
        return (kx,ky,kz),Ek
        ### experimental as not homogenous sampling
        # kx,ky,kz = np.unique(k[:,0]),np.unique(k[:,1]),np.unique(k[:,2])
        # nx,ny,nz = len(kx),len(ky),len(kz)
        # return kx,ky,kz,Ek.reshape((nx,ny,nz,self.nbds))

    def write_Ek(self,Emin=-10,Emax=10):
        
        infile = open(self.path+self.name+'.out','r')
        outfile = open(self.path+self.name+'.Ek','w')
        
        for line in infile:
            line = line.split()
            if line==['Eigenvalues']:
                break
        for _ in range(2):
            next(infile)       
        mu,nbands = self.mu,self.nbds
        nlines = ((np.prod(self.Nk)+1)//2)*(4+self.nbds)
        
        outfile.write('\n');
        for line,_ in zip(infile,range(nlines)):          
            line = line.split()
            if line==[]:
                continue
            elif 'kloop=' in line[0]:
                outfile.write('   '+'   '.join(line)+'\n')
            elif line[::2]==['k1=','k2=','k3=']:
                outfile.write('   '+'   '.join(line)+'\n')
            else:
                Ei = float(line[1])*Ha
                if Ei>(Emin+mu) and Ei<(Emax+mu):
                    outfile.write('   '+'   '.join(line)+'\n')
                if float(line[0])==nbands:
                    outfile.write('   nbands'+'\n\n')        
        infile.close() 

    def bs(self,ext='BANDDAT',spin='1'):
        
        k,E = [],[]; ki,Ei = [],[]
        infile = open(self.path+self.name+'.%s%s'%(ext,spin),'r')    
        for line in infile: 
            line = line.split()
            if line==[]:
                if ki!=[] and Ei!=[]:
                    k.append(ki); E.append(Ei)
                    ki = []; Ei = []
                continue
            else: 
                ki.append(float(line[0]))
                Ei.append(float(line[1]))
        infile.close()
        # k,E 
        for i,ki in enumerate(k):
            if not all(x==ki[0] for x in ki):
                self.Nkpath = i+1
                break
        k,E = k[self.Nkpath-1:],E[self.Nkpath-1:]
        # kpath
        kpts,Nkpts = [],0
        self.kpath = np.zeros(self.Nkpath+1,dtype=int)
        for i in range(self.Nkpath):
            kpts += k[i]; Nkpts += len(k[i])
            self.kpath[i+1] = Nkpts-1
        kpts = np.asarray(kpts)
        # eigenvalues
        Nval = int(len(k)/self.Nkpath)        
        Eval = np.zeros((Nkpts,Nval))
        for i in range(Nval):
            Ei = []
            for j in range(self.Nkpath):
                Ei += E[i*self.Nkpath+j]
            Eval[:,i] = Ei        
        return kpts,Eval

    def DoS(self,path='PDoS'):
        
        file = self.path+path+'/'+self.name+'.DOS.Tetrahedron'
        E,DoS = read(file,num=True),read(file,index=1,num=True)                
        return E,DoS
    
    def PDoS(self,atom='',path='PDoS',name='Tetrahedron',orb=False):
        
        if atom=='':
            atom = range(1,self.Na+1)
        elif isinstance(atom,int) or isinstance(atom,float):
            atom = range(atom,atom+1)
            
        E,PDoS = [],[]
        for i in atom:
            file = self.path+path+'/'+self.name+'.PDOS.%s.atom%s'%(name,i)
            if orb:
                Ei,PDoSi = [],[]
                idx = self.species.index(self.atm[i-1])
                orb_name = [k for k in self.bs[idx][::2]]
                for orb_j in orb_name:
                    if orb_j=='s':
                        j0,j1 = 1,2
                    elif orb_j=='p':
                        j0,j1 = 1,4
                    elif orb_j=='d':
                        j0,j1 = 1,6
                    for j in range(j0,j1):
                        fij = file+'.%s%s'%(orb_j,j)
                        Eij,PDoSij = read(fij,num=True),read(fij,index=1,num=True)
                        Ei.append(Eij); PDoSi.append(PDoSij)
                E.append(Ei); PDoS.append(PDoSi)
            else:
                Ei,PDoSi = read(file,num=True),read(file,index=1,num=True)
                E.append(Ei); PDoS.append(PDoSi)
                
        return np.asarray(E),np.asarray(PDoS) 

    def unfold_orb(self,orb_idx=[],Emin=0,Emax=0,spin='up'):        
        
        if orb_idx==[]:
            orb_idx = np.asarray(self.orb_id).flatten()
        elif not isinstance(orb_idx,np.ndarray):
            orb_idx = np.asarray(orb_idx)
        if Emin==0 and Emax==0:
            Emin,Emax = self.Ebounds
            
        orb_idx += 1 # line = ki Ei orbi0 ...        
        k,E,orb_proj = [],[],[]
        infile = open(self.path+self.name+'.unfold_orb%s'%(spin),'r')    
        for line in infile:
            line = line.split()
            Ei = float(line[1]);
            if Ei > Emin and Ei < Emax :
                k.append(float(line[0]))
                E.append(Ei)
                orbi_proj = [float(line[i]) for i in orb_idx]      
                orb_proj.append(orbi_proj)
        infile.close()      
        return np.asarray(k),np.asarray(E),np.asarray(orb_proj)
        
    def fatbands(self,k_orb,E_orb,orb_proj,dE=.01,DE=.04,Dk=1e-2):
               
        # Intensity map for spectral weight
        # c.f. http://www.openmx-square.org/openmx_man3.9/node171.html
    
        kmin,kmax,nk = min(k_orb),max(k_orb),self.Nk_orb
        # nk = int((kmax-kmin)/dk)+1
        k = np.linspace(kmin,kmax,nk) # kpoints path
        Emin,Emax = self.Ebounds
        nE = int((Emax-Emin)/dE)+1
        E = np.linspace(Emin,Emax,nE) # energy 
        
        _,n_orb = orb_proj.shape
        Imap = np.zeros((nk,nE,n_orb)) # orbital spectral weight 
        X,Y = np.meshgrid(k,E); X,Y = X.T,Y.T
        
        nik,niE = 10*int(round(Dk/(k[1]-k[0]))),10*int(round(DE/dE))
        for i,orbi in enumerate(orb_proj):
            ik,iE = np.argmin(abs(k-k_orb[i])),np.argmin(abs(E-E_orb[i]))
            ik0,ik1,iE0,iE1 = max(0,ik-nik),min(ik+nik,nk),max(0,iE-niE),min(iE+niE,nE)
            Li = 1/(((X[ik0:ik1,iE0:iE1]-k_orb[i])/Dk)**2+((Y[ik0:ik1,iE0:iE1]-E_orb[i])/DE)**2+1)
            for j,orbij in enumerate(orbi):
                Imap[ik0:ik1,iE0:iE1,j] += Li*orbij
        return k,E,Imap
    
    def G0(self,n=0,V=0,SO=False):
        
        id0,id1 = 5,7
        if n==0:
            n = self.Nkt
        if V==0:
            V = abs(np.prod(np.linalg.det(self.acell[1:,1:])))
        if SO:
            id1 = id0
            
        E = self.Et
        if E==np.zeros(0):
            f0 = open(self.path+self.name+'.tran0_0','r')
            for line in (line.split() for line in f0):
                if '#' in line[0]:
                    continue
                else:
                    E.append(float(line[3]))
            E = np.asarray(E)
            self.NEt = len(E)
            
        k = np.zeros((n[0],n[1],2))
        T = np.zeros((n[0],n[1],self.NEt,2))
        for i in range(0,n[0]):
            for j in range(0,n[1]):
                fij = open(self.path+self.name+'.tran%s_%s'%(i,j),'r')
                for line in fij:
                    line = line.split()
                    if line[1::2]==['k2=','k3=']:
                        k[i,j] = np.array([float(line[2]),float(line[4])])
                    elif '#' in line[0]:
                        continue
                    else: 
                        iE = int(line[0])
                        T[i,j,iE,:] = float(line[id0]),float(line[id1])
                fij.close()  
        if SO:
            return k,E,T*.5/V
        return k,E,T/V

###############################################################################
