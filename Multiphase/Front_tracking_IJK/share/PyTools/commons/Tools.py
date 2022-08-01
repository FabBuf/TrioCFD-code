

import matplotlib 
from matplotlib.lines import Line2D
from math import floor, ceil
from matplotlib import rc
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 4}
rc('font', **font)
from TrustFiles import DTEVFile, BuildFromPath
import numpy as np
from cmath import sin
import matplotlib.pyplot as plt
import math
import pickle

#matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['text.latex.unicode'] = True
# matplotlib.rcParams["backend"] = "Qt4Agg"
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{gensymb}",r"\usepackage{txfonts}",r"\usepackage{nicefrac}",r"\usepackage{amssymb}",r"\usepackage{sistyle}"]
#matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{gensymb}",r"\usepackage{txfonts}",r"\usepackage{nicefrac}",r"\usepackage{amssymb}"] #,r"\usepackage{sistyle}"]
#rc('legend', fontsize='medium',numpoints=2)
#rc('text', usetex = True)
#rc('text', usetex = True)
#rc('text', usetex=True)
#rc('font', family='serif')  
# import matplotlib.pyplot as pl
spaceDims = 3
dims = range(spaceDims) 

class Field(object):
    """ Field object is used for scalar, vectorial or tensorial quantity in order to operate between them all
     the operation we want, for instance product of a vector by a matrix or tensorial product between a vector and
     a 3D tensor ans so on."""
    __TRUST_FILES_CACHE = {}
    _VOID_ARR = np.array([], dtype=np.float32)
    #_VOID_MAT2 = np.array([[]], dtype=np.float32)
    global compteur
    compteur=0
    
    def __init__(self, name=None, ax=None, labelRef=None, tex=None):
        """ 
        initialization of a field object
        @param a name for your field
        @param an axis under the form np.array([]) 
        @param a name for this axis (for the plot)
        @param a latex name for the legend
        """
        self.labelRef=labelRef
        self.name = name  # Nom d'usage
        self.tex = tex  # Nom latex
        self._ax = ax  # Axis
        self._npa = np.array([])        
        self.label = np.array([])     
    
    def setname(self, col):
        self.name=col
        
    def settex(self, col):
        self.tex=col
    
    def setlabelref(self, col):
        self.labelRef=col
                
    def __add__(self, other):
        """ 
        @param re-define of the addition operator +
        """
        try:
            res=self.add(other, '+')
        except:
            print 'BE CAREFULL : An + operator has been reversed'
            res=other.add(self, '+')
        return res
    
    def __sub__(self, other):
        """ 
        @param re-define of the addition operator -
        """
        try:
            res=self.add(other, '-')
        except:
            print 'BE CAREFULL : An - operator has been reversed'
            res=other.add(self, '-')
        return res   
        

    def __mul__(self, other):
        """ 
        @param re-define of the product operator *
        """
        #try:
        res=self.product(other, '*')
        #except:
            #print 'BE CAREFULL : An * operator has been reversed'
            #res=other.product(self,'*')
        return res        
        
   
    def createField(self,name):
        """ 
        creation of a Field from an other Field
        @param the Field you want to copy 
        @param a name for the new field
        """
        if isinstance(self,Tensor2Field):
            grad=Tensor2Field(name, self._ax, self.labelRef)
        elif isinstance(self,ScalarField):
            grad=ScalarField(name, self._ax, self.labelRef)
        elif isinstance(self,Tensor3Field):
            grad=Tensor3Field(name, self._ax, self.labelRef)            
        elif isinstance(self,VectorField):
            grad=VectorField(name, self._ax, self.labelRef) 
        else:
            print 'BE CAREFULL : A Field object has been created with no specific scalar/vector or tensor attribute'
            grad=Field(name+self.name, self._ax, self.labelRef) 
        return grad 

    
    def TypeField(self):
        """ 
        Return the type of self
        @param a Field
        """ 
        if isinstance(self,Tensor2Field):
            a='tensor2'
        elif isinstance(self,ScalarField):
            a='scalar'
        elif isinstance(self,Tensor3Field):
            a='tensor3'           
        elif isinstance(self,VectorField):
            a='vector' 
        return a     
    
    
    def compo(self, compo_i, compo_j=None):
        """ 
        Return one component of a second order tensor
        @param a Field (2nd order tensor)
        """  
        a=self.TypeField()
        if a=='vector':
            other=ScalarField()               
            ab = len(self._npa[0,:])
            other._npa=np.zeros((1,ab))
            other.name='bla'
            other._ax=self._ax
            other._npa[0,:]=self._npa[compo_i,:]
            other.labelRef=self.labelRef
            try:
                other.settex(self.tex)
            except:
                pass        
        
        if a=='tensor2':
            other=ScalarField()               
            ab = len(self._npa[compo_i,compo_j,:])
            other._npa=np.zeros((1,ab))
            other.name='bla'
            other._ax=self._ax
            other._npa[0,:]=self._npa[compo_i, compo_j,:]
            other.labelRef=self.labelRef
            try:
                other.settex(self.tex) #+" (%s,%s)"%(compo_i+1,compo_j+1))
            except:
                pass
        return other
    
    def grad(self, bc=0):
        """ 
        Return the gradient of self
        @param a Field
        """    
        other=self.derivate(bc)
        a=self.TypeField()
        if a=='scalar':
            grad=VectorField()            
            aa = len(self._npa[0,:])
            grad.name= 'grad('+self.name+')' 
            grad._npa=np.zeros((3,aa))
            grad._ax=self._ax
            grad._npa[2,:]=other._npa[0,:]
            grad.labelRef=self.labelRef
            
        elif a=='vector':
            other=self.derivate()   
            grad=Tensor2Field(other.name, other.labelRef, other.tex)
            grad.name='grad('+self.name+')'
            c=len(other._ax[0,:])
            grad._npa=np.zeros((3,3,c))
            grad._ax=other._ax
            grad.labelRef=other.labelRef          
            grad._npa[0,2,:]=other._npa[0,:]
            grad._npa[1,2,:]=other._npa[1,:]
            grad._npa[2,2,:]=other._npa[2,:]
            
        elif a=='tensor2':
            grad=Tensor3Field(other.name, other.labelRef, other.tex)
            grad.name='grad('+self.name+')'
            c=len(other._ax[0,:])
            grad._npa=np.zeros((3,3,3,c))
            grad._ax=other._ax
            grad.labelRef=other.labelRef          
            grad._npa[:,:,2 ,:]=other._npa[:,:,:]
            
        elif a=='tensor3':
            grad=Tensor4Field(other.name, other.labelRef, other.tex)
            grad.name='grad('+self.name+')'
            c=len(other._ax[0,:])
            grad._npa=np.zeros((3,3,3,3,c))
            grad._ax=other._ax
            grad.labelRef=other.labelRef          
            grad._npa[:,:,:,2,:]=other._npa[:,:,:,:]            
            
        return grad
    
    
    def integ_const(self, inter=None):
        """ 
        Return the integrated value
        @param a Field
        """        
        
        dx=self._ax[0,2]-self._ax[0,1]
        som=0.0
        
        for i in range(len(self._npa[:,0]) ):
            for j in range(len(self._npa[i,:]) ):
                if j>inter[0] and j<inter[1]:
                    som=som+self._npa[i,j]*dx
        return som    
    
    def integ(self):
        """ 
        Return the integrated field
        @param a Field
        """ 
        inte=self.createField('/['+self.name+']d'+self.labelRef)
        a = len(self._npa[0,:])
        b = len(self._npa[:,0])
        inte._npa=np.zeros((b,a))
        inte._ax=self._ax
        
        for i in range(0, len(self._npa[:,0]) ):
            somme=0
            for j in range(0, len(self._npa[i,:])-1 ):  
                somme=somme+self._npa[i,j]*(self._ax[0,j+1]-self._ax[0,j])
                inte._npa[i,j]=somme
    
        inte._npa[:,-1]=2*inte._npa[:,-2]-inte._npa[:,-3]
        inte.tex=self.tex
        return inte
    
    def tronc(self, val=0.0):
        other=self*0.0
        for i in range(len(self._npa[2,:])):
            if self._npa[2,i]>val:
                other._npa[2,i]=0.0  
            else:
                other._npa[2,i]=self._npa[2,i]
        
        return other   
    
        
    def MoyenneGlissante(self, db=0.4, maxcurv=0.0, wall=0.1):
        """ 
        Return the integrated field
        @param a scalar Field
        """ 
        inte=self.createField('/['+self.name+']d'+self.labelRef)
        a = len(self._npa[0,:])
        b = len(self._npa[:,0])
        inte=self*0.0
        cc=0.0
        for i in range(len(self._npa[0,:])):
            iter=0
            inte._npa[:,i]=0.0
            for j in range(len(self._npa[0,:])-1 ):
                if self._ax[0,j]>wall and  self._ax[0,j]<self._ax[0,i] and self._ax[0,j]<wall+db :
                    cc=cc+1.0
                    if cc==1.0:
                        val=self._npa[0,i]
                    inte._npa[:,i]=inte._npa[:,i]+self._npa[0,j]
                    iter=iter+1.0
                if self._ax[0,i]<=wall*1.3:
                    inte._npa[:,i]=maxcurv
                    
            if iter !=0 : 
                inte._npa[:,i]=inte._npa[:,i]/iter


        inte._npa[:,-1]=inte._npa[:,-2] 
        inte.tex=self.tex
        return inte
        
    def derivate(self, bc_type=0):
        """ 
        Return the derivative of self
        @param a Field
        @param for other boundary conditions
        """ 
        grad=self.createField('d('+self.name+')/d'+self.labelRef)       
        a = len(self._npa[0,:])
        b = len(self._npa[:,0])
        grad._npa=np.zeros((b,a))
        i=0
        while i<=len(self._npa[:,0])-1:
            grad._npa[i,:]= cell_to_cell_gradient_scalar(self._npa[i,:], self._ax[0,:], bc_type)
            i+=1
              
        return grad 
    
    def inv(self, tol=0.0):
        """ 
        Return the inverse of self
        @param a Field
        """ 
        grad=self.createField('inv')      
        a = len(self._npa[0,:])
        b = len(self._npa[:,0])
        grad._npa=np.zeros((b,a))
        for i in range(a):
            for j in range(b):
                if abs(self._npa[j,i])<tol:
                    grad._npa[j,i]=0.0
                    raise Exception('Warning inverse of %s below tolerance %g'%(self.name_, tol))
                else:
                    grad._npa[j,i]=1.0/self._npa[j,i]

        return grad 
    
    def sqrt(self):
        """ 
        Return the square root of self
        @param a Field
        """ 
        grad=self.createField('sqrt')      
        grad._npa=np.sqrt(self._npa)
        for i in range(3):   
            for j in range(3):
                for k in range(len(self._npa[0,0,:]) ):                 
                    if self._npa[i,j,k]>0:
                        grad._npa[i,j,k]=np.sqrt(self._npa[i,j,k])
                    else:
                        grad._npa[i,j,k]=-np.sqrt(-self._npa[i,j,k])
        return grad 

 
 
 
    def symetriser(self, symet=False,axe=None):
        """ 
        Return the symetric of Field 
        @param a field
        @param symet=1 for symetry. -1 for antisymetry. Default is False
        """  
                #### utilise
        
	if isinstance(self,Tensor2Field):
            self.preReshape()

        sym=self*(1.0)
        a = len(self._npa[0,:])
        b = len(self._npa[:,0])
        if axe==None:
            axe=range(b)
        
        
        for i in axe:
            if symet==False:
                sym._npa[i,:]=self._npa[i,:]
            elif symet[i]==1:
                sym._npa[i,:]=(self._npa[i,:]+self._npa[i,::-1])/2.
            elif symet[i]==-1:         
                sym._npa[i,:]=(self._npa[i,:]-self._npa[i,::-1])/2. 
                        
                        
        self._npa=sym._npa
	if isinstance(self,Tensor2Field):
            self.postReshape()

        return self 
    
    

    def contract(self,ref,indice):
        """ 
        allow to do contraction on Field object
        @param a field
        @param ref=0, 'i' or 'ij' for an out-scalar, vector or tensor
        @param the indicies listing for contraction. ex: 'ikjk'
        
        """  
        list=[]
        ind_indep=0
        for i in indice:
            ind_indep=ind_indep+1
            if i in list:
                ind_indep=ind_indep-2
            list.append(i)
        if ref==0:
            pass
        elif ind_indep>len(ref) :
            raise Exception('there is to much independent indicies')
        elif ind_indep<len(ref) :
            raise Exception('there is less independent indicies that expected')            
                
        if ref==0:
            product=ScalarField(self.name, self._ax, self.labelRef)
            c=len(self._ax[0,:])
            product._npa=np.zeros((1,c))        
            for k in dims:
                if indice=='kk':
                    product._npa[0,:]=product._npa[0,:]+self._npa[k,k,:]
                elif indice=='klkl':
                    for l in dims:                   
                        product._npa[0,:]=product._npa[0,:]+self._npa[k,l,k,l,:]
		else:
		    raise Exception("Contraction case %s -> %s not coded"%(ref,indice)  )   
            
              
        elif len(ref)==1:
            product=VectorField(self.name, self._ax, self.labelRef)
            c=len(self._ax[0,:])
            product._npa=np.zeros((3,c))  
            
            for i in range(0, len(self._npa[:,0,0,0]) ):        
                for k in range(0, len(self._npa[0,0,:,0]) ):
                    if indice=='ikk':
                        product._npa[i,:]=product._npa[i,:]+self._npa[i,k,k,:]
		    else:
		    	raise Exception("Contraction case %s -> %s not coded"%(ref,indice)  )  
        elif len(ref)==2:    
            product=Tensor2Field(self.name, self._ax, self.labelRef)
            c=len(self._ax[0,:])
            product._npa=np.zeros((3,3,c))  
            
            for i in range(0, len(self._npa[:,0,0,0]) ):        
                for j in range(0, len(self._npa[0,:,0,0]) ): 
                    for k in range(0, len(self._npa[0,0,:,0]) ):
                        if indice=='kijk':
                            product._npa[i,j,:]=product._npa[i,j,:]+self._npa[k,i,j,k,:]
                        elif indice=='kjik':   
                            product._npa[i,j,:]=product._npa[i,j,:]+self._npa[k,j,i,k,:]
                        elif indice=='ikjk':   
                            product._npa[i,j,:]=product._npa[i,j,:]+self._npa[i,k,j,k,:] 
                        elif indice=='ijkk':   
                            product._npa[i,j,:]=product._npa[i,j,:]+self._npa[i,j,k,k,:]
                        elif indice=='ikkj':   
                            product._npa[i,j,:]=product._npa[i,j,:]+self._npa[i,k,k,j,:]  
                        elif indice=='kikj':   
                            product._npa[i,j,:]=product._npa[i,j,:]+self._npa[k,i,k,j,:]
			else:
			    raise Exception("Contraction case %s -> %s not coded"%(ref,indice)  )  
        elif len(ref)==3: 
            #resultat un Tensor3
            product=Tensor3Field(self.name, self._ax, self.labelRef)
            c=len(self._ax[0,:])
            product._npa=np.zeros((3,3,3,c))  
            for i in range(0, len(self._npa[:,0,0,0]) ):        
                for j in range(0, len(self._npa[0,:,0,0]) ): 
                    for k in range(0, len(self._npa[0,0,:,0]) ): 
                        if indice=='kji': 
                            product._npa[i,j,k,:]=product._npa[i,j,k,:]+self._npa[k,j,i,:] 
                        elif indice=='jik': 
                            product._npa[i,j,k,:]=product._npa[i,j,k,:]+self._npa[j,i,k,:]   
			else:
			    raise Exception("Contraction case %s -> %s not coded"%(ref,indice)  )    
        else:
	    raise Exception("Contraction case %s -> %s not coded"%(ref,indice)  )                                                                          
        return product
        
    
    
    def product(self, other, sym='*'):

        """ 
        Isomorphist product operator 
        @param a field
        @param a field
        """        
        if isinstance(other,float):
            product=self.createField(self.name+sym+str(other))
            product._npa=other*self._npa
            product._ax=self._ax
 

        elif isinstance(other,ScalarField) and isinstance(self,Tensor2Field):
            print 'attention modif recente'
            product=self*0.0
            for i in range(len(product._npa[:,0])):	   
            	product._npa[i,:]=self._npa[i,:]*other._npa[0,:]
        
        else:

            product=other.createField(self.name+sym+other.name)
            aa = len(self._npa[0,:])    
            ab = len(self._npa[:,0])
            da = len(other._npa[0,:])    
            db = len(other._npa[:,0]) 
            product._npa=np.zeros((max(ab,db),max(aa,da)))
            product._ax=self._ax
               
            if db%ab==0 and aa==da:
                dim=db/ab
                for i in range(0, len(other._npa[:,0]) ):
                    a=i/dim
                    for j in range(0, len(other._npa[0,:]) ):
                        if sym=='*':
                            product._npa[i,j]=self._npa[a,j]*other._npa[i,j]
                        if sym=='/':  
                            product._npa[i,j]=other._npa[i,j]/self._npa[a,j]
            else :
		#import pdb; pdb.set_trace()
                raise Exception('Ce produit est impossible a effectuer. Vector lengths mismatch.')
                     
        return product
    
    def add(self,other,sym=None):
        """ 
        Isomorphist add operator 
        @param a field
        @param a field
        """ 
        try:
            som=self.createField('('+self.name+sym+other.name+')')        
        except:
            som=self.createField('('+self.name+sym+str(other)+')')
        if isinstance(other,float):
            if sym=='+':
                som._npa=other+self._npa
            if sym=='-':
                som._npa=self._npa-other
            if sym=='--':
                som._npa=other-self._npa      
        else:
            aa = len(self._npa[0,:])
            ad = len(other._npa[0,:])
            ba = len(self._npa[:,0])
            bd = len(other._npa[:,0])
            som._npa=np.zeros((ba,aa))
            if aa!=ad or ba!=bd:
                raise Exception('ces champs ne peuvent pas etre additiones %d, %d, %d, %d'%(aa, ad, ba, bd))
                return
            
            for i in range(0, len(self._npa[:,0]) ):
                for j in range(0, len(self._npa[0,:]) ):
                    if sym=='+':
                        som._npa[i,j]=self._npa[i,j]+other._npa[i,j]
                    if sym=='-':
                        som._npa[i,j]=self._npa[i,j]-other._npa[i,j]    
        return som

    def ChangerMaillage(self, other):
	tmp=other*0.0
	for i in range(len(tmp._ax[0,:])):
		for j in range(len(self._ax[0,:])):
			if self._ax[0,j]>tmp._ax[0,i]:
				avant=tmp._ax[0,i]-self._ax[0,j-1]
				apres=self._ax[0,j]-tmp._ax[0,i]
				av=j-1
				ap=j
				break
		print "verif maille", avant+apres, av, ap
		tmp._npa[:,i]=(avant*self._npa[:,ap]+apres*self._npa[:,av])/(avant+apres)
	return tmp

    @classmethod
    def LoadFromFile(cls, filename, col_ref, cols, name, labelRef, tex, bc=0):
        #### utilise
        verb=0
        if verb:
           print 'Load %s from file %s waiting .....................'%(name, filename)
        fileObj = cls.getObject(filename)
        yoyo=[]
        x, tmp = fileObj.getValues(col_ref[0])
        #print col_ref[0]
        #print x
        #print tmp
        yoyo.append(tmp)
        
        if len(cols) == 1:
            l=[]
            ret = ScalarField(name,np.vstack(yoyo),labelRef,tex)
            x, tmp = fileObj.getValues(cols[0])
            l.append(tmp)
	    if (bc==2):
            	l[0]=(l[1]+l[-2])*0.5
            	l[-1]=(l[1]+l[-2])*0.5
            ret._npa_x = x
            ret._npa = np.vstack(l)
            ret._label=cols
            
        elif len(cols) == spaceDims:
            ret = VectorField(name,np.vstack(yoyo),labelRef,tex)
            for i in dims:
                x, tmp = fileObj.getValues(cols[i])
                if i == 0:
                    nz = tmp.shape[0]
                    l = np.empty((spaceDims, nz), dtype=float)  
                    l[i,:] = tmp  
	            if (bc==2):
            		l[i,0]=(l[i,1]+l[i,-2])*0.5
            		l[i,-1]=(l[i,1]+l[i,-2])*0.5
                else:
                    if tmp.shape[0] != nz:
                        raise Exception
                    l[i,:] = tmp 
	            if (bc==2):
            		l[i,0]=(l[i,1]+l[i,-2])*0.5
            		l[i,-1]=(l[i,1]+l[i,-2])*0.5
                    pass
            ret._npa_x = x
            ret._npa = np.vstack(l)
            ret._label=cols         
            
            
        elif len(cols) == spaceDims*spaceDims:
            ret = Tensor2Field(name,np.vstack(yoyo),labelRef,tex)
            for i in dims:
                for j in dims:
                    x, tmp = fileObj.getValues(cols[i+spaceDims*j])
                    if i==0 and j==0:
                        nz = tmp.shape[0]
                        l = np.empty((spaceDims, spaceDims, nz), dtype=float)  
                        l[i,j,:] = tmp 
			if (bc==2):
            			l[i,j,0]=(l[i,j,1]+l[i,j,-2])*0.5
            			l[i,j,-1]=(l[i,j,1]+l[i,j,-2])*0.5  
                    else:
                        if tmp.shape[0] != nz:
                            raise Exception
                        l[i,j,:] = tmp 
			if (bc==2):
            			l[i,j,0]=(l[i,j,1]+l[i,j,-2])*0.5
            			l[i,j,-1]=(l[i,j,1]+l[i,j,-2])*0.5 

                    pass
                pass
            
            ret._npa_x = x
            ret._npa=l  
            ret._label=cols 
            
        elif len(cols) == spaceDims*spaceDims*spaceDims:
            ret = Tensor3Field(name,np.vstack(yoyo),labelRef,tex)
            for i in dims:
                for j in dims:
                    for k in dims:
                        x, tmp = fileObj.getValues(cols[i+spaceDims*j+spaceDims*spaceDims*k])
                        if i==0 and j==0 and k==0:
                            nz = tmp.shape[0]
                            l = np.empty((spaceDims, spaceDims, spaceDims, nz), dtype=float)  
                            l[i,j,k,:] = tmp  
			    if (bc==2):
            			l[i,j,k,0]=(l[i,j,k,1]+l[i,j,k,-2])*0.5
            			l[i,j,k,-1]=(l[i,j,k,1]+l[i,j,k,-2])*0.5 
                        else:
                            if tmp.shape[0] != nz:
                                raise Exception
                            l[i,j,k,:] = tmp 
			    if (bc==2):
            			l[i,j,k,0]=(l[i,j,k,1]+l[i,j,k,-2])*0.5
            			l[i,j,k,-1]=(l[i,j,k,1]+l[i,j,k,-2])*0.5 

                    pass
                pass
            
            ret._npa_x = x
            ret._npa=l  
            ret._label=cols     
        
        else:
            raise Exception("Number of columns %d given not recognized" % (len(cols)))
            raise NotImplementedError
            pass
        
        if verb:
           print 'end.'
        return ret
    
    
    @classmethod
    def initgravity(cls, const, name, other):
        """ 
        Return a vector Field initialize with three constance
        @param a list of three scalar number
        @param a name
        @param a 'model Field' for x axis and labels
        """ 
        var=VectorField(name,other._ax,other.labelRef,other.tex)
        c=len(other._ax[0,:])
        var._npa=np.zeros((3,c))
        var._npa[0,:]=const[0]
        var._npa[1,:]=const[1]
        var._npa[2,:]=const[2]
        return var  
          
    @classmethod
    def initsource(cls, const, name, other):      
        """ 
        Return a scalar field
        @param scalar number
        @param a name
        @param a 'model Field' for x axis and labels
        """ 
        var=VectorField(name,other._ax,other.labelRef,other.tex)
        c=len(other._ax[0,:])
        var._npa=np.zeros((3,c))
        var._npa[0,:]=const[0]
        var._npa[1,:]=const[1]
        var._npa[2,:]=const[2]
        return var  
       
    @classmethod   
    def initFromTab(cls, genre,u, v=0, w=0):
        """ a suprimer """
        if genre=='vec':
            var=VectorField()
            var.setname('bla')
            var._ax=np.zeros((1,len(u[:,0])))
            var._npa=np.zeros((3,len(u[:,0])))
            var._ax[0,:]=u[:,0]
            var._npa[0,:]=u[:,1]
            var._npa[1,:]=v[:,1]
            var._npa[2,:]=w[:,1]
        if genre=='sca':
            var=ScalarField()
            var.setname('bla')
            var._ax=np.zeros((1,len(u[:,0])))
            var._npa=np.zeros((1,len(u[:,0])))
            var._ax[0,:]=u[:,0]
            var._npa[0,:]=u[:,1]
        return var 
     
    def initFromAnalytic(cls, genre,a, b, N, name, label, labelRef, tex):
        """ 
        initialize field with analytical expression
        @param genre could be S, V, T, T3 or T4 for scalar, vector, Tensor, 3rd order tensor and 4th order tensor
        @param a=x(0)
        @param b=x(N)
        @param N=number of point
        @param a name corresponding to a register case
        @param a label
        @param a label for the x axis
        @param a name for the legend
        """ 
        if genre=='S':
            var=ScalarField(name,np.linspace(a,b,N).reshape(1,N),labelRef,tex)
            c=len(var._ax[0,:])
            if name=='zero':
                var._npa=0*np.sin(var._ax)
            if name=='un':
                var._npa=np.sin(var._ax)/np.sin(var._ax)                 
            if name=='sin':
                var._npa=np.sin(var._ax)
            if name=='2sin':
                var._npa=2*np.sin(var._ax)  
            if name=='sin**2':
                var._npa=np.sin(var._ax)*np.sin(var._ax)                   
            if name=='cos':
                var._npa=np.cos(var._ax)    
            if name=='tan':
                var._npa=np.tan(var._ax)             
        if genre=='V':
            var=VectorField(name,np.linspace(a,b,N).reshape(1,N),labelRef,tex)
            c=len(var._ax[0,:])
            var._npa=np.zeros((3,c))
            if name== 'poiseuille':
                var._npa[0,:]=-0.5*(2*var._ax-var._ax*var._ax)
            if name=='poly':   
                var._npa[0,:]=var._ax[:]
                var._npa[1,:]=0.5*var._ax*var._ax
                var._npa[2,:]=0.166666*var._ax*var._ax*var._ax
            if name=='dpoly':   
                var._npa[0,:]=1
                var._npa[1,:]=var._ax[:]
                var._npa[2,:]=0.5*var._ax*var._ax         
            if name=='titi':
                var._npa[0,:]=1
                var._npa[1,:]=np.cos(var._ax)
                var._npa[2,:]=np.sqrt(np.abs(var._ax))
        if genre=='T':
            var=Tensor2Field(name,np.linspace(a,b,N).reshape(1,N),labelRef,tex)
            c=len(var._ax[0,:])
            var._npa=np.zeros((3,3,c))
            if name=='identite':
                a=np.zeros((3,c))
                a[0,0,:]=1.0
                a[1,1,:]=1.0
                a[2,2,:]=1.0   
            if name=='VV':
                a=np.zeros((3,c))
                a[0,:]=var._ax[:]
                a[1,:]=0.5*var._ax*var._ax
                a[2,:]=0.166666*var._ax*var._ax*var._ax
                for i in dims:
                    for j in dims:
                        var._npa[i,j,:]=a[i,:]*a[j,:]
            if name=='sin':
                for i in dims:
                    for j in dims:
                        var._npa[i,j,:]=np.sin(var._ax)
            if name=='const':
                var._npa[0,0,:]=1
                var._npa[1,0,:]=2
                var._npa[2,0,:]=3
                var._npa[0,1,:]=4
                var._npa[1,1,:]=5
                var._npa[2,1,:]=6
                var._npa[0,2,:]=7
                var._npa[1,2,:]=8
                var._npa[2,2,:]=9
            if name== 'transpo':
                var._npa[0,0,:]=1
                var._npa[0,1,:]=2
                var._npa[0,2,:]=3
                var._npa[1,0,:]=4
                var._npa[1,1,:]=5
                var._npa[1,2,:]=6
                var._npa[2,0,:]=7
                var._npa[2,1,:]=8
                var._npa[2,2,:]=9
            if name=='dconst':
                var._npa[0,0,:]=0
                var._npa[1,0,:]=0
                var._npa[2,0,:]=0
                var._npa[0,1,:]=0
                var._npa[1,1,:]=0
                var._npa[2,1,:]=0
                var._npa[0,2,:]=0
                var._npa[1,2,:]=0
                var._npa[2,2,:]=0
            if name== 'poiseuille':
                var._npa[0,2,:]=-0.5*(2-2*var._ax)
                        
            if name== 'ijkk':
                for i in dims:
                    for j in dims:
                        var._npa[i,j,:]=(i+1)+j*3+(i+1)+j*3+9+27+(i+1)+j*3+18+54
            if name== 'kijk':
                for i in dims:
                    for j in dims:
                        var._npa[i,j,:]=1+i*3+j*9+2+i*3+j*9+27+3+i*3+j*9+2*27
            if name== 'ikjk':
                for i in dims:
                    for j in dims:
                        var._npa[i,j,:]=1+j*3+i*9+2+j*3+i*9+27+3+j*3+i*9+2*27
            if name== 'ikjk':
                for i in dims:
                    for j in dims:
                        var._npa[i,j,:]=(i+1)+j*9+(i+1)+3+j*9+27+(i+1)+2*3+j*9+2*27                                               

        if genre=='T3':
            var=Tensor3Field(name,np.linspace(a,b,N).reshape(1,N),labelRef,tex)
            c=len(var._ax[0,:])
            var._npa=np.zeros((3,3,3,c))
            if name=='gradsin':
                for i in dims:
                    for j in dims:
                        var._npa[i,j,2,:]=np.cos(var._ax)
            if name=='sin':
                for i in dims:
                    for j in dims:
                        for k in dims:
                            var._npa[i,j,k,:]=np.sin(var._ax)  
                            
            if name=='TV':
                a=np.zeros((3,c))
                a[0,:]=var._ax[:]
                a[1,:]=0.5*var._ax*var._ax
                a[2,:]=0.166666*var._ax*var._ax*var._ax
                for i in dims:
                    for j in dims:
                        for k in dims:
                            for l in dims:
                                var._npa[i,j,k,:]=((i+1)+j*3)*a[k,:]
            if name=='VT':
                a=np.zeros((3,c))
                a[0,:]=var._ax[:]
                a[1,:]=0.5*var._ax*var._ax
                a[2,:]=0.166666*var._ax*var._ax*var._ax
                for i in dims:
                    for j in dims:
                        for k in dims:
                            for l in dims:
                                var._npa[i,j,k,:]=((j+1)+k*3)*a[i,:]                               
                                
            if name=='file':
                for i in dims:
                    for j in dims:
                        for k in dims:
                            var._npa[i,j,k,:]=(i+1)+j*3+k*9                                                
                           
        if genre=='T4':
            var=Field(name,np.linspace(a,b,N).reshape(1,N),labelRef,tex)
            c=len(var._ax[0,:])
            var._npa=np.zeros((3,3,3,3,c))
            if name=='cte':
                for i in dims:
                    for j in dims:
                        for k in dims:
                            for l in dims:
                                var._npa[i,j,k,l,:]=(i+1)+j*3+k*9+l*27
                
            if name=='gradsin':
                for i in dims:
                    for j in dims:
                        for k in dims:
                            var._npa[i,j,k,2,:]=np.cos(var._ax)
            if name=='TT':
                for i in dims:
                    for j in dims:
                        for k in dims:
                            for l in dims:
                                var._npa[i,j,k,l,:]=((i+1)+j*3)*((k+1)+l*3)
            
            if name=='T3V':
                a=np.zeros((3,c))
                a[0,:]=var._ax[:]
                a[1,:]=0.5*var._ax*var._ax
                a[2,:]=0.166666*var._ax*var._ax*var._ax
                for i in dims:
                    for j in dims:
                        for k in dims:
                            for l in dims:
                                var._npa[i,j,k,l,:]=np.sin(var._ax)*a[l,:]                                           
        return var

    @classmethod
    def getObject(cls, filename):
        fileObj = cls.__TRUST_FILES_CACHE.get(filename, None)
        if fileObj is None:
            fileObj = DTEVFile(filename, None)
            cls.__TRUST_FILES_CACHE[filename] = fileObj
            pass
        return fileObj

    @classmethod
    def getEntries(cls, filename):
        fObj = cls.getObject(filename)
        dvar = fObj.getEntries()
        #print "Available variables are : ", dvar
        return dvar

    @classmethod
    def printEntries(cls, filename):
        dvar = cls.getEntries(cls, filename)
        print "Available variables are : ", dvar
        return dvar

        

class ScalarField(Field):
    def __init__(self, name=None, ax=None, labelRef=None, tex=None):
        Field.__init__(self, name, ax,labelRef, tex)
        
    def product(self, other, sym='*'):

        """ 
        return kind of matrix product. But be carefull with the definition.
        Usefull for product of a tensor with a scalar. For complex operation use 
        TensorProduct() and contract()
        @param a Field 
        @param a Field
        """  

        if isinstance(other,Tensor2Field) and other!=self:
            other.preReshape()
            resultat = Field.product(self,other,sym)
        else:
            resultat = Field.product(self,other,sym)
            
        if isinstance(other,Tensor2Field) and other!=self:
            other.postReshape()
        if isinstance(resultat,Tensor2Field):
            resultat.postReshape()
        return resultat
    
    def mean(self):
        """ 
        Return the mean value of self
        @param a Field (scalar Field)
        """         
        resultat=np.mean(self._npa[:])
        return resultat  
    
    def sqrt(self):
        """ 
        Return the square root of self
        @param a Field (scalar Field)
        """ 
        grad=self.createField('sqrt')      
        grad._npa=np.sqrt(self._npa)  
        for j in range(len(self._npa[0,:])):                 
            if self._npa[0,j]>0:
                grad._npa[0,j]=np.sqrt(self._npa[0,j])
            else:
                grad._npa[0,j]=-np.sqrt(-self._npa[0,j])
        return grad 
    
    def power(self, p=2.0):
        """ 
        Return self to the power p
        @param a Field (scalar Field)
        @param a scalar number
        """ 
        grad=self*0.0    
        for j in range(len(self._npa[0,:])): 
            #print self._npa[0,j], p, self._npa[0,j]**p
            grad._npa[0,j]=self._npa[0,j]**p
        return grad  

    def tanh(self):
        """ 
        Return tanh(self)
        @param a Field (scalar Field)
        @param a scalar number
        """ 
        grad=self*0.0    
        for j in range(len(self._npa[0,:])): 
            #print self._npa[0,j], p, self._npa[0,j]**p
            grad._npa[0,j]=np.tanh(self._npa[0,j])
        return grad 

    def exp(self):
        """ 
        Return self to the power p
        @param a Field (scalar Field)
        @param a scalar number
        """ 
        grad=self*0.0    
        for j in range(len(self._npa[0,:])): 
            grad._npa[0,j]=np.exp(self._npa[0,j])
        return grad  

    
    def abs(self):
        """ 
        Return self to the power p
        @param a Field (scalar Field)
        @param a scalar number
        """ 
        grad=self*0.0    
        for j in range(len(self._npa[0,:])): 
            #print self._npa[0,j], p, self._npa[0,j]**p
            grad._npa[0,j]=abs(self._npa[0,j])
        return grad      
     
    
    def ValueInX(self, x):  
        dx=self._ax[0,2]-self._ax[0,1]
        for i in range(len(self._npa[0,:])):
            if self._ax[0,i]>x-dx and self._ax[0,i]<x+dx:
                val= self._npa[0,i]
        return val

    def IndxInX(self, x):  
        dx=self._ax[0,2]-self._ax[0,1]
        for i in range(len(self._npa[0,:])):
            if self._ax[0,i]>x-dx and self._ax[0,i]<x+dx:
                a=i
        return a   
    
    def ProjectionChampsDifferentesTaille(self, other):
        #on veut projeter les valeurs de self sur l axe de other
        cut=ScalarField()               
        ab = len(other._npa[0,:])
        cut._npa=np.zeros((1,ab))
        cut.name='bla'
        cut._ax=other._ax
        cut.labelRef=self.labelRef
        try:
            cut.settex(self.tex)
        except:
            pass        
        for i in range(len(cut._npa[0,:])):
            val=-1
            for j in range(len(self._ax[0,:])):
                if self._ax[0,j]<cut._ax[0,i] and self._ax[0,j+1]>cut._ax[0,i]:
                    val=j
            if val==-1:
                print 'il a pas trouver de correspondance entre les deux maillages'
            else:
                cut._npa[0,i]=self._npa[0,val] 
        return cut
    
    
class VectorField(Field):
    def __init__(self, name=None, ax=None, labelRef=None, tex=None):
        Field.__init__(self, name, ax,labelRef, tex)


    def mean(self):
        """ 
        Return the mean value of self
        @param a Field (scalar Field)
        """         
	resultat=self*0.0
        resultat._npa[0,:]=np.mean(self._npa[0,:])
        resultat._npa[1,:]=np.mean(self._npa[1,:])
        resultat._npa[2,:]=np.mean(self._npa[2,:])
	return resultat

    def mean_compo(self, compo):
        """ 
        Return the mean value of self
        @param a Field (scalar Field)
        """         
	return np.mean(self._npa[compo,:])
        
    def MatProduct(self,other):
        """ a supprimer a terme par TensorProduct"""
        b=other.TypeField()
        if b =='tensor2':
            pro=self.product(other)
            res=self.createField(self.name+'*'+other.name)            
            aa = len(self._npa[0,:])    
            ab = len(self._npa[:,0])
            res._npa=np.zeros((ab,aa))
            res._npa[0,:]=pro._npa[0,0,:]+pro._npa[0,1,:]+pro._npa[0,2,:]
            res._npa[1,:]=pro._npa[1,0,:]+pro._npa[1,1,:]+pro._npa[1,2,:]    
            res._npa[2,:]=pro._npa[2,0,:]+pro._npa[2,1,:]+pro._npa[2,2,:]
            
        elif  b=='vector': 
            res=Tensor2Field()            
            aa = len(self._npa[0,:])    
            ab = len(self._npa[:,0])
            res._npa=np.zeros((ab,ab,aa))
            res._ax=self._ax
            res.name=self.name+'.'+other.name
            res.labelRef='z'
            for i in dims:
                for j in dims:
                    res._npa[i,j,:]=self._npa[i,:]*other._npa[j,:]
          
        return res
    

    def product(self, other, sym='*'):

        """ 
        return kind of matrix product. But be carefull with the definition.
        Usefull for product of a tensor with a scalar. For complex operation use 
        TensorProduct() and contract()
        @param a Field 
        @param a Field
        """    
        if isinstance(other,Tensor2Field) and other!=self:
            other.preReshape()
            resultat = Field.product(self,other,sym)
            other.postReshape()
            resultat.postReshape()
        elif isinstance(other,Tensor3Field) and other!=self:
            resultat = Field.vecMatProduct(self,other,sym)
        elif isinstance(other,VectorField):
            resultat = Field.product(other,self,sym)
        elif isinstance(other,ScalarField):
            resultat = Field.product(other,self,sym)
        else:
            resultat = Field.product(self,other,sym)
                  
        return resultat

    def power(self, p=2.0):
        """ 
        Return self to the power p
        @param a Field (scalar Field)
        @param a scalar number
        """ 
        grad=self*0.0    
        for i in range(3):
        	for j in range(len(self._npa[i,:])): 
            	#print self._npa[0,j], p, self._npa[0,j]**p
            		grad._npa[i,j]=self._npa[i,j]**p
        return grad 
 
    def sqrt(self):
        """ 
        Return self to the power p
        @param a Field (scalar Field)
        @param a scalar number
        """ 
        grad=self*0.0    
        for i in range(3):
        	for j in range(len(self._npa[i,:])): 
            	#print self._npa[0,j], p, self._npa[0,j]**p
                        if self._npa[i,j]>0:
            			grad._npa[i,j]=np.sqrt(self._npa[i,j])
			else:
            			grad._npa[i,j]=-np.sqrt(-self._npa[i,j])			
        return grad  
    
    def tensorProduct(self, other):
        """ 
        return tensorial product between a vectir and a Field.
        @param a Field (1rd order tensor)
        @param a Field (Vector or Tensor2 Field)
        """      

        c=len(other._ax[0,:])
        if isinstance(other, Tensor2Field):
            product=Tensor3Field(self.name+'.'+str(other))
            product._npa=np.zeros((3,3,3,c))
        if isinstance(other, VectorField):
            product=Tensor2Field(self.name+'.'+str(other))
            product._npa=np.zeros((3,3,c))
                
        product._ax=self._ax
        product.labelRef=other.labelRef 
        for i in dims:   
            for j in dims:
                if isinstance(other, VectorField):
                    product._npa[i,j,:]=self._npa[i,:]*other._npa[j,:]
                if isinstance(other, Tensor2Field):
                    for k in dims:
                        product._npa[i,j,k,:]=self._npa[i,:]*other._npa[j,k,:]
        
        return product


class Tensor2Field(Field):
    def __init__(self, name=None, ax=None, labelRef=None, tex=None):
        Field.__init__(self, name, ax,labelRef, tex)
    
    def preReshape(self):
        """ 
        Take a tabular 3*3*n and give a 9*n  
        @param a Field 
        """  
        a = len(self._npa[0,0,:])
        self._npa=self._npa.reshape(9,a) 
        
    def postReshape(self):
        """ 
        Take a tabular 9*n and give a 3*3*n  
        @param a Field 
        """  
        a = len(self._npa[0,:])
        self._npa=self._npa.reshape(3,3,a)
        
    def derivate(self, bc_type=0):
        """ 
        Give the derivative  
        @param a Field 
        @param bc_type for other boundary conditions
        
        """  
        self.preReshape()
        resultat = Field.derivate(self, bc_type)
        self.postReshape()
        resultat.postReshape()
        return resultat
    
   
    
    def add(self,other,sym):
        """ 
        return self+other
        @param a Field 
        @param a Field
        """  
        try :
            self.preReshape()
            if other!=self:
                other.preReshape()
        except :
            print 'addition avec un scalaire'    
        resultat = Field.add(self,other,sym)
        self.postReshape()
        if isinstance(other, float):
            pass
        elif other!=self:
            other.postReshape()
        resultat.postReshape()
        return resultat
    
    def product(self, other, sym='*'):

        """ 
        return kind of matrix product. But be carefull with the definition.
        Usefull for product of a tensor with a scalar. For complex operation use 
        TensorProduct() and contract()
        @param a Field 
        @param a Field
        """          


        if isinstance(other,Tensor2Field) and other!=self:
            self.preReshape()
            other.preReshape()
            resultat = Field.product(self,other,sym)
        if isinstance(other,VectorField):
            self.preReshape()
            resultat = Field.product(other,self,sym)
        if isinstance(other,float):
            resultat = Field.product(self,other,sym)
        if isinstance(other,ScalarField):
            self.preReshape()          
            resultat = Field.product(self,other,sym)

        if isinstance(other,Tensor2Field) and other!=self:
            self.postReshape()
            other.postReshape()
        if isinstance(other,VectorField):
            self.postReshape()
        if isinstance(other,ScalarField):
            self.postReshape()  
        if isinstance(resultat,Tensor2Field) and isinstance(other,float)==0:
            resultat.postReshape()
        return resultat
    
    def MatProduct(self,other):
        """ a supprimer a terme par TensorProduct"""
        b=other.TypeField()
        if  b=='tensor2':
            #produit mat mat donne mat 
            res=self.createField(self.name+'.'+other.name)            
            aa = len(self._npa[0,:,0])    
            ab = len(self._npa[:,0,0])
            ac = len(self._npa[0,0,:])
            res._npa=np.zeros((ab,aa,ac))
            for k in range(len(res._npa[0,0,:])):
                res._npa[:,:,k]=np.dot(self._npa[:,:,k],other._npa[:,:,k])
        return res
    
    def transpo(self):
        """ 
        return the transpose tensor
        @param a Field (2rd order tensor)
        """      
        
        res=self.createField(self.name+'^Tr')
        aa = len(self._npa[0,:,0])    
        ab = len(self._npa[:,0,0])
        ac = len(self._npa[0,0,:])
        res._npa=np.zeros((ab,aa,ac))
        for k in range(len(self._npa[0,0,:])):  
            res._npa[:,:,k]=np.transpose(self._npa[:,:,k])
        return res
    
    
    def tensorProduct(self, other):
        """ 
        return tensorial product between a second order tensor and a Field (Vector orTensor2).
        @param a Field (2rd order tensor)
        @param a Field (Vector or Tensor Field)
        """      

        c=len(other._ax[0,:])
        
        if isinstance(other, Tensor2Field):
            product=Tensor4Field(self.name+'.'+str(other))
            product._npa=np.zeros((3,3,3,3,c))
        elif isinstance(other, VectorField):
            product=Tensor3Field(self.name+'.'+str(other))
            product._npa=np.zeros((3,3,3,c)) 
        else:
	    raise Exception("tensorProduct other type %s not coded"%(other.type())  )                
               
        #product=self.createField(self.name+'.'+str(other))
        product._ax=self._ax
        product.labelRef=other.labelRef 
        
        for i in dims:   
            for j in dims:
                for k in dims:
                    if isinstance(other, VectorField):
                        product._npa[i,j,k,:]=self._npa[i,j,:]*other._npa[k,:]
                    if isinstance(other, Tensor2Field):
                        for l in dims:
                            product._npa[i,j,k,l,:]=self._npa[i,j,:]*other._npa[k,l,:]
    
                           
        return product

    def tensoradd(self, other):
        """ 
        return the add between two second order tensor.
        @param a Field (2rd order tensor)
        @param a Field (2rd order tensor)
        """        
        c=len(other._ax[0,:])
        
        if isinstance(other, Tensor2Field):
            product=Tensor2Field(self.name+'.'+str(other))
            product._npa=np.zeros((3,3,c))

        product._ax=self._ax
        product.labelRef=other.labelRef 
        
        for i in dims:   
            for j in dims:
                for k in dims:
                    product._npa[i,j,:]=product._npa[i,j,:]+self._npa[i,k,:]+other._npa[k,j,:]
    
                           
        return product
    
class Tensor3Field(Field):
    def __init__(self, name=None, ax=None, labelRef=None, tex=None):
        Field.__init__(self, name, ax,labelRef, tex)
        
    def preReshape(self):
        """ 
        Take a tabular 3*3*3*n and give a 27*n  
        @param a Field 
        """  
        a = len(self._npa[0,0,0,:])
        self._npa=self._npa.reshape(27,a) 
        
    def postReshape(self):
        """ 
        Take a tabular 27*n and give a 3*3*3*n  
        @param a Field 
        """  
        a = len(self._npa[0,:])
        self._npa=self._npa.reshape(3,3,3,a)
        
    def derivate(self, bc_type=0):
        """ 
        Give the derivative  
        @param a Field 
        @param bc_type for other boundary conditions
        
        """  
        self.preReshape()
        resultat = Field.derivate(self, bc_type)
        self.postReshape()
        resultat.postReshape()
        return resultat
    
    def add(self,other,sym):
        """ 
        return self+other
        @param a Field 
        @param a Field
        """  
        try :
            self.preReshape()
            if other!=self:
                other.preReshape()
        except :
            print 'addition avec un scalaire'    
        resultat = Field.add(self,other,sym)
        self.postReshape()
        if isinstance(other, float):
            pass
        elif other!=self:
            other.postReshape()
        resultat.postReshape()
        return resultat
    
    def product(self, other, sym='*'):

        """ 
        return kind of matrix product. But be carefull with the definition.
        Usefull for product of a tensor with a scalar. For complex operation use 
        TensorProduct() and contract()
        @param a Field 
        @param a Field
        """  
        self.preReshape()
        if isinstance(other,Tensor2Field) and other!=self:
            other.preReshape()
            resultat = Field.product(self,other,sym)
        if isinstance(other,VectorField):
            resultat = Field.product(other,self,sym)
        if isinstance(other,float):
            resultat = Field.product(self,other,sym)
        
        self.postReshape()
        if isinstance(other,Tensor2Field) and other!=self:
            other.postReshape()
        if isinstance(resultat,Tensor2Field):
            resultat.postReshape()
        if isinstance(resultat,Tensor3Field):
            resultat.postReshape()   
        return resultat
    
#     def MatProduct(self,other):      
#         b=other.TypeField()
#         if  b=='tensor2':
#             #produit mat mat donne mat 
#             res=self.createField(self.name+'.'+other.name)            
#             aa = len(self._npa[0,:,0])    
#             ab = len(self._npa[:,0,0])
#             ac = len(self._npa[0,0,:])
#             res._npa=np.zeros((ab,aa,ac))
#             for k in range(len(res._npa[0,0,:])):
#                 res._npa[:,:,k]=np.dot(self._npa[:,:,k],other._npa[:,:,k])
#         return res
    # GB ajout du scalar to be checked
    def tensorProduct(self, other):
        """ 
        return tensorial product between a third order tensor and a vector.
        @param a Field (3rd order tensor)
        @param a Field (Scalar or Vector Field)
        """  

        c=len(other._ax[0,:])
        if isinstance(other, VectorField):
            product=Tensor4Field(self.name+'.'+str(other))
            product._npa=np.zeros((3,3,3,3,c))
	elif isinstance(other, ScalarField):
            product=Tensor3Field(str(other)+self.name)
            product._npa=np.zeros((3,3,3,c))
	else:
	    raise Exception("Error in %s::tensorProduct"%self.TypeField())

        product._ax=self._ax
        product.labelRef=other.labelRef 
        
        for i in dims:   
            for j in dims:
                for k in dims:
                    if isinstance(other, ScalarField):
		        product._npa[i,j,k,:]=self._npa[i,j,k,:]*other._npa[:]
                    elif isinstance(other, VectorField):
                        for l in dims:
                            product._npa[i,j,k,l,:]=self._npa[i,j,k,:]*other._npa[l,:]
                           
        return product

class Tensor4Field(Field):
    def __init__(self, name=None, ax=None, labelRef=None, tex=None):
        Field.__init__(self, name, ax,labelRef, tex)
    
    pass
#############################""

def fluctuij(p, du, pdu):
    """ 
    return the mean field of the correlation fluctuation of two variables.
    @param mean of A
    @param mean of B
    @param mean of AB
    return mean of A'B'
    """    
    if isinstance(p, VectorField) and isinstance(du, VectorField):
        pdu_bis=p.MatProduct(du)
    if isinstance(p, ScalarField) and isinstance(du, VectorField):
        pdu_bis=p*du
    if isinstance(p, ScalarField) and isinstance(du, Tensor2Field):   
        pdu_bis=du*p
    if isinstance(p, Tensor2Field) and isinstance(du, Tensor2Field):    
        pdu_bis=p.tensorProduct(du)
        pdu_bis=pdu_bis.contract('ij', 'ikjk')
    res=pdu-pdu_bis
    return res

def fluctuijk(uuu, u, uu):
    """ 
    return the mean field of the correlation fluctuation of three variables.
    @param mean of AAA
    @param mean of A
    @param mean of AA
    return mean of A'A'A'
    """    
    Rijk=uuu
    Rij_uk=u.tensorProduct(uu)
    #print Rij_uk
    uijk=u.tensorProduct(u).tensorProduct(u)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                Rijk._npa[i,j,k,:]=Rijk._npa[i,j,k,:]+uijk._npa[i,j,k,:]*2.0-Rij_uk._npa[j,i,k,:]-Rij_uk._npa[i,j,k,:]-Rij_uk._npa[k,i,j,:] 
         
    #Rijk._npa[0,0,0,:]=0.0
    return Rijk    

def fluctuij_with_I(p, du, pdu,I):

    """ 
    return the mean field of the correlation fluctuation of two variables.
    @param mean of A
    @param mean of B
    @param mean of AB
    return mean of A'B'
    """    
    if isinstance(p, VectorField) and isinstance(du, VectorField):
        pdu_bis=p.MatProduct(du)
        pdu_bis=pdu_bis*(I.inv(0.0001))
    if isinstance(p, ScalarField) and isinstance(du, VectorField):
        pdu_bis=p*du*(I.inv(0.0001))

    if isinstance(p, ScalarField) and isinstance(du, Tensor2Field):   
        pdu_bis=du*p*(I.inv(0.0001))

    if isinstance(p, Tensor2Field) and isinstance(du, Tensor2Field):    
        pdu_bis=p.tensorProduct(du)
        pdu_bis=pdu_bis.contract('ij', 'ikjk')
        pdu_bis=pdu_bis*(I.inv(0.0001))
    res=pdu-pdu_bis
    return res

def fluctuijk_with_I(uuu, u, uu, I):

    """ 
    return the mean field of the correlation fluctuation of three variables.
    @param mean of AAA
    @param mean of A
    @param mean of AA
    return mean of A'A'A'`

    """    
    Rijk=uuu
    Rij_uk_bis=u.tensorProduct(uu) # GB: I think it is in fact ui_Rjk ?

    #print Rij_uk
    uijk_bis=u.tensorProduct(u).tensorProduct(u)
    if isinstance(I, ScalarField) and isinstance(Rij_uk_bis, Tensor3Field) and isinstance(uijk_bis, Tensor3Field): 
        Rij_uk=Rij_uk_bis*1.
        uijk = uijk_bis*1.

    for i in range(3):
        for j in range(3):
            for k in range(3):
                # GB : For a true parmutation, it should not be Rij_uk._npa[j,i,k,:] but Rij_uk._npa[j,k,i,:] 
                # But it doesn't matter because the permutation (i,k) -> (k,i) plays on the Rij part (last 2 indices) that are symmetric.
                Rijk._npa[i,j,k,:]=Rijk._npa[i,j,k,:]+uijk._npa[i,j,k,:]*2.0-Rij_uk._npa[j,i,k,:]-Rij_uk._npa[i,j,k,:]-Rij_uk._npa[k,i,j,:] 
         
    #Rijk._npa[0,0,0,:]=0.0
    return Rijk    

def load(fichier):
    """ 
    load a dictionnary  
    @param a file name
    """  
    output=open(fichier, 'rb')
    Run=pickle.load(output)
    output.close()
    clef(Run)
    return Run

def loadSansClef(fichier):
    """ 
    load a dictionnary without keys  
    @param a file name
    """  
    output=open(fichier, 'rb')
    Run=pickle.load(output)
    output.close()
    return Run

def save(Var, fichier):
    """ 
    save a dictionnary  
    @param a dictionnary  
    @param a file name
    """       
    
    output=open(fichier, 'wb')
    pickle.dump(Var,output)
    output.close()
    return

def clef(SphericalRun):
    """ 
    print all the keys of a dictionnary  
    @param a dictionnary  
    """   
    print ''
    print 'Liste des objets'
    for cle in SphericalRun.keys():
        i=0
        try:
            for cle2 in SphericalRun[cle].keys():
                i=i+1
                if i==1:
                    print cle+'--->'+cle2
                else:
                    print '------->'+cle2 
        except:
                print '--->'+cle
           
    print ''
    return


def NetContribution(Mk, db=0.3, L=2.0):
    """ 
    return mean value without taking into account 0. 
    @param the Field (vector field)
    """
    dk=Mk._ax[0,2]-Mk._ax[0,1]
    MkNet=Mk*0.0
    for j in range(len(Mk._npa[:,0])):
        for i in range(len(Mk._npa[j,:])):
            som=0.0
            somax=0.0
            #### trouver les bornes 
            Ln=len(Mk._npa[j,:])
            dbn=int(Ln*db/L)  
            indix=[]
            for k in range(int(dbn/2)):
                if k==0 :
                    indix.append(i-k)
                    continue
                if i-k in range(Ln):
                    indix.append(i-k)
                if i+k in range(Ln):
                    indix.append(i+k)
            
            
            indix.sort()
            for k in indix[:-1]:
                    som=som+(0.5*Mk._npa[j,k]+0.5*Mk._npa[j,k+1])
                    somax=somax+1.0
                    

            #print i, len(indix), len(indix[:-1]), somax, som, Mk._npa[j,i], som+Mk._npa[j,i]
            if somax==0:
                MkNet._npa[j,i]=0.0
            else:
                MkNet._npa[j,i]=som/somax
            #print i, somax
    return MkNet

def absoluteValue(Mk):
    """ 
    return the absolute value
    @param the Field (vector field)
    """
    MkNet=Mk*0.0
    for j in range(len(Mk._npa[:,0])):
        for i in range(len(Mk._npa[j,:])):
            MkNet._npa[j,i]=abs(Mk._npa[j,i])
    return MkNet

def cell_to_cell_gradient_scalar(Field, dField, bc_type=0):
    """ 
    return the first derivative of Field 
    @param the Field. Component by component
    @param the axis field
    @param bc_type for other boundary condition
    """
    gradi = np.zeros(len(Field))
  
    gradi[1:-1] = (Field[2:]-Field[:-2])/(dField[2:]-dField[:-2])
    if bc_type==1:
        for i in range(len(gradi[1:-1])/4):
            if i==0:
                gradi[i]=(-0.166666666*Field[i+6]+1.2*Field[i+5]-3.75*Field[i+4]+6.66666666*Field[i+3]-7.5*Field[i+2]+6*Field[i+1]-2.45*Field[i])/(dField[i+1]-dField[i])
            else:
                gradi[i]=(-0.166666666*Field[i+6]+1.2*Field[i+5]-3.75*Field[i+4]+6.66666666*Field[i+3]-7.5*Field[i+2]+6*Field[i+1]-2.45*Field[i])/(dField[i+1]-dField[i])
                gradi[-i]=(0.166666666*Field[-i-6]-1.2*Field[-i-5]+3.75*Field[-i-4]-6.66666666*Field[-i-3]+7.5*Field[-i-2]-6*Field[-i-1]+2.45*Field[-i])/(dField[-i]-dField[-i-1])
       
    elif bc_type==0:
        gradi[0]=(-0.166666666*Field[6]+1.2*Field[5]-3.75*Field[4]+6.66666666*Field[3]-7.5*Field[2]+6*Field[1]-2.45*Field[0])/(dField[1]-dField[0])
        gradi[-1]=(0.166666666*Field[-7]-1.2*Field[-6]+3.75*Field[-5]-6.66666666*Field[-4]+7.5*Field[-3]-6*Field[-2]+2.45*Field[-1])/(dField[-1]-dField[-2])
	
    elif bc_type==2: #condition periodique
        gradi[0] = (Field[1]-Field[-1])/(dField[1]-dField[-1])
        gradi[-1] = (Field[0]-Field[-2])/(dField[0]-dField[-2])

    return gradi   

def FromListToScalarField(liste, name="bla", labelRef="x", tex="x"):
	# Before : 
	var=ScalarField(name,(np.array(liste[0])).reshape(1,len(liste[0])),labelRef,tex)
	var._npa=np.array(liste[1]).reshape(1,len(liste[1]))
	# GB modif 2020.04.01 : 
	#var=ScalarField(name,(np.array(liste[0])).reshape(1,len(liste[0]))[0],labelRef,tex)
	#var._npa=np.array(liste[1]).reshape(1,len(liste[1]))[0]
	return var

def autoscale_y(ax,margin=0.1):
    """This function rescales the y-axis based on the data that is visible given the current xlim of the axis.
    ax -- a matplotlib axes object
    margin -- the fraction of the total height of the y-data to pad the upper and lower ylims"""
    def get_bottom_top(line):
        xd = line.get_xdata()
        yd = line.get_ydata()
        lo,hi = ax.get_xlim()
        y_displayed = yd[((xd>lo) & (xd<hi))]
        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed)-margin*h
        top = np.max(y_displayed)+margin*h
        return bot,top

    lines = ax.get_lines()
    bot,top = np.inf, -np.inf

    for line in lines:
        new_bot, new_top = get_bottom_top(line)
        if new_bot < bot: bot = new_bot
        if new_top > top: top = new_top

    ax.set_ylim(bot,top)
    return ax


class FictiveBubble(object):
    def __init__(self,amplitude=0.0, istart=0, Ist=0.0, dIst=0.0, ifin=0, If=0.0, dIf=0.0,ipic=0, Ip=0.0, dIp=0.0, alpha=0.0, fax=0.0, iax=0.0, pax=0.0 ):
        self.Start=[istart, Ist, dIst, iax]
        self.flatStart=[istart, Ist, dIst, iax]
        self.flatEnd=[istart, Ist, dIst, iax]
        self.fin=[ifin, If, dIf, fax]
        self.pic=[ipic, Ip, dIp, pax]
        self.alpha=alpha
        self._npa = np.array([])        
        self.label = np.array([])
        self.amplitude=amplitude
        
    def findPic(self, I):
        maxi=max(abs(self._npa[:]-1))
        for i in range(len(self._npa[:])):
            if abs(self._npa[i]-1)==maxi:
                self.pic[0]=i+self.Start[0]
                self.pic[1]=self._npa[i]
                self.pic[2]=I.grad()._npa[0,i+self.Start[0]]
                self.pic[3]=I._ax[0, i+self.Start[0]]
                
        print 'pic de taux de vide =', 1-self.pic[1], 'soit x=', self.pic[3], 'et dI/dz=', self.pic[2]
    
        
