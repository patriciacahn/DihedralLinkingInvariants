from dihedrallinking import *



print('Knot: 10_121')
p=5
print('p=',p)
q=int((p+1)/2)% p

print('The local writhe numbers at the heads of arcs k0,k1,...,k9 are:')
signlist=[-1, 1, -1, -1, -1, 1, -1, -1, 1, -1]
print(signlist)

print('The subscripts of the arcs passing over the heads of arcs k0,k1,...,k9 are:')
overstrands=[3, 8, 9, 6, 2, 1, 0, 4, 5, 7]
print(overstrands)

print('The colors of the arcs k0,k1,...,k9 are:')
colorlist=[1, 2, 3, 4, 3, 3, 1, 1, 5, 1]
print(colorlist)


arrowlists=createarrowlists(p,colorlist,overstrands)
vertexlists=createvertexlists(p,colorlist,overstrands)
for j in range(q):
    for k in range(q):
        print('Intersection number of K^j, j=', j, 'with Sigma^k, k=',k )
        print(intKjSigmak(j, k, p, overstrands, signlist, arrowlists, vertexlists))