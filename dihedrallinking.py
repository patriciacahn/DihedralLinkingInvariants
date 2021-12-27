
def reflect(i,j,p):
    #reflect vertex i over line through j on p gon
    s=(2*j-i)%p
    return s


def createarrowlists(myp,mycolorlist,myoverstrands):
    n=len(mycolorlist)
    myq=int((myp+1)/2) %myp
    vertices=[]
    v=mycolorlist[0]
    for i in range(myq):
        vertices.append(v)
        v=(v+1)%myp    
    myinitiallist=[]
    for i  in vertices:
        myinitiallist.append([reflect(i,mycolorlist[0],myp),i])
    myarrowlists=[]
    myarrowlists.append(myinitiallist)
    for i in range(n-1):              
        newlist=[]
        for j in range(myq):           
            newlist.append([reflect(myarrowlists[i][j][0],mycolorlist[myoverstrands[i]],myp),reflect(myarrowlists[i][j][1],mycolorlist[myoverstrands[i]],myp)])
        myarrowlists.append(newlist)    
    return  myarrowlists


def createvertexlists(myp,mycolorlist,myoverstrands):

    
    n=len(mycolorlist)    
    myvertexlists=[]
    myarrowlists=createarrowlists(myp,mycolorlist,myoverstrands)
    for i  in range(n):
        myvertexlists.append(arrowtovertexlist(myp,myarrowlists[i]))
    return myvertexlists


def arrowtovertexlist(myp,myarrowlist) :
    myq=int((myp+1)/2) %myp
    myvertexlist=[ [ 0 for i in range(2) ] for j in range(myp) ]
    for i in range(myq) :  
        tail=myarrowlist[i][0]
        head=myarrowlist[i][1]
        if head==tail:
            myvertexlist[head][0]=i
            myvertexlist[head][1]='c'
        else:
            myvertexlist[tail][0]=i
            myvertexlist[tail][1]='t'
            myvertexlist[head][0]=i
            myvertexlist[head][1]='h'
    return myvertexlist

def above(i,j,myp,myoverstrands,myarrowlists,myvertexlists):
    s= myvertexlists[myoverstrands[i]][myarrowlists[i][j][1]][0]
    return s

def below(i,j,myp,myoverstrands,myarrowlists,myvertexlists):
    s= myvertexlists[myoverstrands[i]][myarrowlists[i][j][0]][0]
    return s

def epsilona(i,j,myp,myoverstrands,myarrowlists,myvertexlists):
    a=above(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
    if myarrowlists[i][j][1]==myarrowlists[myoverstrands[i]][a][1] and a!=0:       
        return 1
    elif myarrowlists[i][j][1]==myarrowlists[myoverstrands[i]][a][0] and a!=0:
        return -1
    elif a==0:
        return 0
    
def epsilonb(i,j,myp,myoverstrands,myarrowlists,myvertexlists):
    b=below(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
    if myarrowlists[i][j][0]==myarrowlists[myoverstrands[i]][b][0] and b!=0:
        return 1
    elif myarrowlists[i][j][0]==myarrowlists[myoverstrands[i]][b][1] and b!=0:
        return -1
    elif b==0:
        return 0
    
def Ca(i,j,k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists) :
    a=above(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
    epsa=epsilona(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
    if a==k and k!=0 and mysignlist[i]*epsa==-1:
        return -mysignlist[i]
    elif a==k and k==0:
        return -mysignlist[i]
    else:
        return 0
    
def Cb(i,j,k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists) :
    b=below(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
    epsb=epsilonb(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
    if b==k and k!=0 and mysignlist[i]*epsb==1:
        return mysignlist[i]
    elif b==k and k==0:
        return mysignlist[i]
    else:
        return 0 
    
def matrix2chain(k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists):
    n=len(myoverstrands)
    myq=int((myp+1)/2) %myp
    dim=(myq-1)*n
    M=[[0 for i in range(dim+1)] for j in range(dim)]
    for i in range(n):
        for j in range(1,myq):
            a=above(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
            b=below(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
            epsa=epsilona(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
            epsb=epsilonb(i,j,myp,myoverstrands,myarrowlists,myvertexlists)
            ca=Ca(i,j,k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists)
            cb=Cb(i,j,k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists)
            M[(myq-1)*i+j-1][n*(j-1)+i]+=1
            M[(myq-1)*i+j-1][n*(j-1)+((i+1)%n)]+=-1
            if a!=0:
                M[(myq-1)*i+j-1][n*(a-1)+myoverstrands[i]]+=-epsa
            if b!=0:
                M[(myq-1)*i+j-1][n*(b-1)+myoverstrands[i]]+=-epsb
            M[(myq-1)*i+j-1][dim]+=-ca-cb
            
    
            
    return M

def coeflist(k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists):
    n=len(myoverstrands)
    myq=int((myp+1)/2) %myp
    dim=(myq-1)*n
    M=matrix2chain(k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists)
    M=Matrix(M)
    R=M.rref()
    pivots=list(R[1])
    R=R[0]
    constants=R.col(dim).tolist()
    coefs=[0 for x in range(dim)]
    for x in range(len(pivots)):
      coefs[pivots[x]]=constants[x][0]
    return coefs
from sympy import *

def intKjSigmak(j,k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists):
    n=len(mysignlist)
    coefs=coeflist(k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists)
    lk=0
    for i in range(n):
        a=above(i, j, myp, myoverstrands, myarrowlists, myvertexlists)
        if a!=0:           
            lk+=epsilona(i, j, myp, myoverstrands, myarrowlists, myvertexlists)*coefs[(a-1)*n+myoverstrands[i]]
        lk+=-Ca(i,j,k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists)
    return lk

def intersectionlist(i,j,k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists):
    n=len(mysignlist)
    coefs=coeflist(k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists)
    a=above(i, j, myp, myoverstrands, myarrowlists, myvertexlists)
    
    intersection=0
    if a!=0:
        intersection+=epsilona(i, j, myp, myoverstrands, myarrowlists, myvertexlists)*coefs[(a-1)*n+myoverstrands[i]]
    intersection+=-Ca(i,j,k,myp,myoverstrands,mysignlist, myarrowlists,myvertexlists)
    return intersection


