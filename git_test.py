import numpy as np
import time
import math
import pickle
from itertools import product

import sage.all
from sage.rings.integer_ring import ZZ
from sage.homology import simplicial_complex as sc





#simple_root = [(0,-1,1,0,0,0,0,0),(0,0,-1,1,0,0,0,0),(0,0,0,-1,1,0,0,0),(0,0,0,0,-1,1,0,0),(0,0,0,0,0,-1,1,0),(0,0,0,0,0,0,-1,1),(1/2,1/2,1/2,1/2,-1/2,-1/2,-1/2,-1/2)]
#fundamental_coweight = [((1,-1,0,0,0,0,0,0),(2,-1,-1,0,0,0,0,0),(3,-1,-1,-1,0,0,0,0),(3,0,0,0,0,1,1,1),(2,0,0,0,0,0,1,1),(1,0,0,0,0,0,0,1),(3,1,1,1,1,1,1,1))]
#name = 'E7'
#dim = 7




#simple_root = [(1,-1,0,0,0,0,0,0),(0,1,-1,0,0,0,0,0),(0,0,1,-1,0,0,0,0),(0,0,0,1,-1,0,0,0),(0,0,0,0,1,-1,0,0),(0,0,0,0,0,1,1,0),(-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2),(0,0,0,0,0,1,-1,0)]
#fundamental_coweight =[((1,0,0,0,0,0,0,-1),(1,1,0,0,0,0,0,-2),(1,1,1,0,0,0,0,-3),(1,1,1,1,0,0,0,-4),(1,1,1,1,1,0,0,-5),(1,1,1,1,1,1,1,-7),(0,0,0,0,0,0,0,-1),(1,1,1,1,1,1,-1,-5))]
#name = 'E8'
#dim = 8



#simple_root_sub = simple_root[1:]
#group = coset_representation(simple_root, fundamental_coweight)










def reflection(u,v) : # return 'the integral primitive vector' that the direction is same with the vector which is reflected 'v' by the hyperplane perpendicular to 'u'
    a = np.dot(u,v)
    b = np.dot(u,u)
    n = len(u)
    T = []
    for i in range(n) :
        T.append(v[i]-2*(a/b)*u[i])
    c = 2
    L = []
    for i in range(n) :
        L.append(int(T[i]))
    while not T == L :
        T = [] 
        for i in range(n) :
            T.append(c*v[i]-2*c*(a/b)*u[i])
        c += 1
        L = []
        for i in range(n) :
            L.append(int(T[i]))
    d = math.gcd(int(T[0]),int(T[1]))
    for i in range(2,n) :
        d = math.gcd(d,int(T[i]))
    S = []
    for i in range(n) :
        S.append(int((1/d)*int(T[i])))
    return S





#coset representation

def coset_representation(simple_root, fundamental_coweight):
    T = []
    group = ['']
    S = [fundamental_coweight[0][0]]
    while T != T+S :
        T = T+S
        basket = []
        groupyb = []
        for i in range(len(simple_root)) :
            for j in range(len(S)) :
                basket.append(tuple(reflection(simple_root[i],S[j])))
                groupyb.append('{}'.format(group[T.index(S[j])])+'{}'.format(i))
        S = []
        for k in range(len(basket)) :
            if basket[k] not in T and basket[k] not in S :
                S.append(basket[k])
                group.append(groupyb[k])
    return group





# Find the vertices of the Coxeter complex from the fundamental coweights
def vtx(simple_root,fundamental_coweight) :
    my_list = []
    for i in range(len(fundamental_coweight[0])) :
        T = []
        S = [fundamental_coweight[0][i]]
        while not T == T+S :
            T = T+S
            basket = set()
            for i in range(len(simple_root)) :
                for j in range(len(S)) :
                    basket.add(tuple(reflection(simple_root[i],S[j])))
            S = list(basket-set(T))
        my_list = my_list + T
    return my_list




def vtxpre(simple_root,fundamental_coweight) :
    my_list = []
    for i in range(len(fundamental_coweight[0])) :
        T = []
        S = [fundamental_coweight[0][i]]
        while not T == T+S :
            T = T+S
            basket = set()
            for i in range(len(simple_root)) :
                for j in range(len(S)) :
                    basket.add(tuple(reflection(simple_root[i],S[j])))
            S = list(basket-set(T))
        my_list.append(T)
    return my_list

def find_index(my_set,vector) :
    my_index = my_set.index(vector)+1
    return my_index

def find_vector(my_set,number) :
    my_vector = my_set[number-1]
    return my_vector





# Find the coordinate of the vector in the coweight lattice
def char_coordinate(vector,simple_root) :
    S = []
    for i in range(len(simple_root)) :
        a = np.dot(simple_root[i],vector)
        S.append(a)
    c = 2
    L = []
    for i in range(len(simple_root)) :
        L.append(int(S[i]))
    while not S == L :
        S = [] 
        for i in range(len(simple_root)) :
            a = np.dot(simple_root[i],vector)
            S.append(c*a)
        c += 1
        L = []
        for i in range(len(simple_root)) :
            L.append(int(S[i]))
    d = math.gcd(S[0],S[1])
    for i in range(2,len(simple_root)) :
        d = math.gcd(d,int(S[i]))
    T = []
    for i in range(len(simple_root)) :
        T.append((int((1/d)*S[i]))%2)
    return T





def find_reduced_rowspace(simple_root) :
    mtx = []
    for i in range(len(simple_root)) :
        S = []
        for j in range(len(simple_root)) :
            S.append(int((2*(np.dot(simple_root[i],simple_root[j]))/(np.dot(simple_root[i],simple_root[i])))%2))
        mtx.append(S)
    S = []
    for i in product([1,0], repeat=len(simple_root)):
        a = i
        S.append(a)
    zero = []
    for i in range(len(simple_root)) :
        zero.append(0)
    S.remove(tuple(zero))
    my_list = []
    while not S == [] :
        T = set()
        set_1 = {S[0]}
        while not T == T.union(set_1) :
            set_2 = set()
            T = T.union(set_1)
            for i in range(len(mtx)) :
                for x in set_1 :
                    if np.dot(mtx[i],x)%2 == 1 :
                        P = list(x)
                        P[i] = (P[i]+1)%2
                        set_2.add(tuple(P))
            set_1 = set_2 - T
        my_list.append([S[0],len(T),T])
        for x in T :
            S.remove(x)
    return my_list




def find_simplex(simple_root,fundamental_coweight) :
    T = set()
    S = list(fundamental_coweight)
    while not T == T.union(set(S)) :
        T = T.union(set(S))
        spx = set()
        for i in range(len(simple_root)) :
            for x in S :
                P = []
                for j in range(len(fundamental_coweight[0])) :
                    P.append(tuple(reflection(simple_root[i],x[j])))
                spx.add(tuple(P))
        S = list(spx - T)
    return list(T)




# return the orbit types and their orbit size
def Type(simple_root):
    reduced_rowspace = find_reduced_rowspace(simple_root)
    Type = []
    T = []
    Num = []
    for i in reduced_rowspace:
        T.append(i[0])
        Num.append(i[1])
    Type.append(T)
    Type.append(Num)
    return Type








# Basic version


def maximal_simplex_basic(simple_root, fundamental_coweight):
    facet = find_simplex(simple_root, fundamental_coweight)
    vtxset = vtx(simple_root,fundamental_coweight)
    row_space = find_reduced_rowspace(simple_root)
    for j in range(len(row_space)) :
        x = row_space[j][0]
        T = set()
        for y in facet :
            S = []
            for i in range(len(y)) :
                S.append((np.dot(x,char_coordinate(y[i],simple_root))%2)*(vtxset.index(y[i])+1))
            while 0 in S:
                         S.remove(0)
            T.add(tuple(S))
        T = list(T)
        T.sort()
        count = 1
        while count != 0 :
            count = 0
            S = []
            for i in range(len(T)-1) :
                if not set(T[i]).issubset(set(T[i+1])) :
                    S.append(T[i])
                else :
                    count +=1
            S.append(T[len(T)-1])
            P = set(S)
            P = list(P)
            P.sort()
            T = P
        P= []
        for i in range(dim+1) :
            P.append([])
        for x in T :
            c = len(x)
            P[c-1].append(x)
        S = set()
        i = 0
        while len(P[i+1]) != 0 :
            for m in range(len(P[i])) :
                x = P[i][m]
                for k in range(len(P[i+1])) :
                    y = P[i+1][k]
                    if set(x).issubset(set(y)) :
                        S.add(x)
                        break
            i += 1
        lst = list(set(T)-S)
        T.sort()
        with open('simplices of {} of type {} basic ver.pkl'.format(name,j+1), 'wb') as h:
            pickle.dump(lst, h, pickle.HIGHEST_PROTOCOL)
        print('Type_{} complete'.format(j+1))





# Find the maximal simplex list using by coset representation


def maximal_simplex_coset(simple_root, fundamental_coweight):
    start = time.time()
    simple_root_sub = simple_root[1:]
    sim_sub = find_simplex(simple_root_sub,fundamental_coweight)
    group = coset_representation(simple_root, fundamental_coweight)
    Typ = Type(simple_root)
    As = []
    Bs = []
    Cs = []
    Ds = []
    for ty in Typ[0] :
        vtxsetp = vtx(simple_root_sub,fundamental_coweight)
        A = []
        B = []
        C = []
        D = []
        for i in range(len(group)) :
            ele = group[i]
            y = fundamental_coweight[0][0]
            for m in range(len(ele)) :
                y = tuple(reflection(simple_root[int(ele[m])],y))
            if (np.dot(char_coordinate(y,simple_root),ty))%2 == 1 :
                A.append(i)
                B.append({i})
                break
        for i in range(len(group)) :
            ele = group[i]
            y = fundamental_coweight[0][0]
            for m in range(len(ele)) :
                y = tuple(reflection(simple_root[int(ele[m])],y))
            if (np.dot(char_coordinate(y,simple_root),ty))%2 == 0 :
                C.append(i)
                D.append({i})
                break
        if len(A)!=0:
            g = group[A[0]]
        if len(C)!=0:
            h = group[C[0]]
        for i in range(len(group)) :
            ele = group[i]
            y = fundamental_coweight[0][0]
            for m in range(len(ele)) :
                y = tuple(reflection(simple_root[int(ele[m])],y))
            if (np.dot(char_coordinate(y,simple_root),ty))%2 == 1 :
                count = 0
                for x in vtxsetp :
                    z = tuple(list(x))
                    for m in range(len(ele)) :
                        z = tuple(reflection(simple_root[int(ele[m])],z))
                    for m in range(len(g)) :
                        x = tuple(reflection(simple_root[int(g[m])],x))
                    if (np.dot(char_coordinate(x,simple_root),ty))%2 == (np.dot(char_coordinate(z,simple_root),ty))%2 :
                        count +=1
                if not count==len(vtxsetp) :
                    A.append(i)
                    g = group[A[-1]]
                    B.append({i})
                else :
                    B[-1].add(i)
            elif (np.dot(char_coordinate(y,simple_root),ty))%2 == 0 :
                count = 0
                for x in vtxsetp :
                    z = tuple(list(x))
                    for m in range(len(ele)) :
                        z = tuple(reflection(simple_root[int(ele[m])],z))
                    for m in range(len(h)) :
                        x = tuple(reflection(simple_root[int(h[m])],x))
                    if (np.dot(char_coordinate(x,simple_root),ty))%2 == (np.dot(char_coordinate(z,simple_root),ty))%2 :
                        count +=1
                if not count==len(vtxsetp) :
                    C.append(i)
                    h = group[C[-1]]
                    D.append({i})
                else :
                    D[-1].add(i)
        As.append(A)
        Bs.append(B)
        Cs.append(C)
        Ds.append(D)
    print(As,Bs,Cs,Ds)
    print("Natural isomorphism : complete")
    typ = Typ[0]
    facet = sim_sub
    for num in range(1,len(typ)+1) :
        start_sub = time.time()
        ty = typ[num-1]
        L = [[As[num-1],Bs[num-1]],[Cs[num-1],Ds[num-1]]]
        for N in L:
            A = N[0]
            B = N[1]
            for numb in range(len(A)) :
                i = A[numb]
                T = set()
                ele = group[i]
                for k in range(len(facet)) :
                    x = facet[k]
                    P = []
                    for j in range(len(fundamental_coweight[0])) :
                        y = x[j]
                        for m in range(len(ele)) :
                            y = reflection(simple_root[int(ele[m])],y)
                        P.append(tuple((np.dot(ty,char_coordinate(y,simple_root))%2)*y))
                    while () in P :
                        P.remove(())
                    T.add(tuple(P))
                lst = list(T)
                S = []
                P = []
                for j in range(dim) :
                    P.append([])
                for j in range(len(lst)) :
                    if len(lst[j]) > 1 :
                        P[len(lst[j])-2].append(lst[j])
                par = 0
                for j in range(len(P)) :
                    par += len(P[j])**2
                co = -1
                while not par == 0 :
                    co += 1
                    par = par - len(P[co])**2
                for j in range(co) :
                    for k in range(len(P[j])) :
                        count = 0
                        for m in range(len(P[j+1])):
                            if set(P[j][k]).issubset(set(P[j+1][m])) :
                                break
                            else :
                                count +=1
                        if count == len(P[j+1]) :
                            S.append(P[j][k])
                lst = list(set(S).union(set(P[co])))
                with open('simplices of {} of type_{} No.{} small.pkl'.format(name,num,i),'wb') as h:
                    pickle.dump(lst, h, pickle.HIGHEST_PROTOCOL)
                B[numb].remove(i)
                for numbe in B[numb] :
                    elep = group[numbe]
                    c = 0
                    if not len(ele) == 0 :
                        while elep[c] == ele[c] and c < len(ele)-1 and c < len(elep)-1 : 
                            c += 1
                    my_str = []
                    for j in range(len(ele)-c) :
                        my_str.append(ele[len(ele)-j-1])
                    for j in range(c,len(elep)) :
                        my_str.append(elep[j])
                    my_lst = []
                    for j in range(len(lst)) :
                        E = []
                        for k in range(len(lst[j])) :
                            z = lst[j][k]
                            for m in range(len(my_str)) :
                                z = tuple(reflection(simple_root[int(my_str[m])],z))
                            E.append(z)
                        my_lst.append(tuple(E))
                    with open('simplices of {} of type_{} No.{} small.pkl'.format(name,num,numbe), 'wb') as h:
                        pickle.dump(my_lst, h, pickle.HIGHEST_PROTOCOL)
        LST = set()
        for j in range(len(group)):
            with open('simplices of {} of type_{} No.{} small.pkl'.format(name,num,j),'rb') as h:
                lst = pickle.load(h)
            LST = LST.union(set(lst))
        with open('simplices of {} of type_{} small.pkl'.format(name,num), 'wb') as h:
            pickle.dump(list(LST), h, pickle.HIGHEST_PROTOCOL)
        print('type_{} complete'.format(num), 'Size of list :', len(LST))




# type_number : int
def number_of_components(simple_root, fundamental_coweight, type_number):
    simple_root_sub = simple_root[1:]
    ty = Type(simple_root)[0][type_number-1]
    vtxsetp = vtx(simple_root_sub, fundamental_coweight)
    group = coset_representation(simple_root, fundamental_coweight)
    cone_label = []
    cone = []
    for i in range(len(group)) :
        ele = group[i]
        y = fundamental_coweight[0][0]
        for m in ele :
            y = tuple(reflection(simple_root[int(m)],y))
        if (np.dot(char_coordinate(y,simple_root),ty))%2 == 1 :
            cone_label.append(i)
            cone.append({i})
            break
    if len(cone_label)!=0:
            g = group[cone_label[0]]
    for i in range(len(group)) :
        ele = group[i]
        y = fundamental_coweight[0][0]
        for m in range(len(ele)) :
            y = tuple(reflection(simple_root[int(ele[m])],y))
        if (np.dot(char_coordinate(y,simple_root),ty))%2 == 1 :
            count = 0
            for x in vtxsetp :
                z = tuple(list(x))
                for m in range(len(ele)) :
                    z = tuple(reflection(simple_root[int(ele[m])],z))
                for m in range(len(g)) :
                    x = tuple(reflection(simple_root[int(g[m])],x))
                if (np.dot(char_coordinate(x,simple_root),ty))%2 == (np.dot(char_coordinate(z,simple_root),ty))%2 :
                    count +=1
            if not count==len(vtxsetp) :
                cone_label.append(i)
                g = group[cone_label[-1]]
                cone.append({i})
            else :
                cone[-1].add(i)
    A = []
    C = []
    for i in cone:
        C += list(i)
    if len(C) == 0:
        print("there is no cone")
        return 1
    for i in C:
        ele = group[i]
        D = []
        for j in range(len(vtxsetp)):
            y = vtxsetp[j]
            for m in ele:
                y = tuple(reflection(simple_root[int(m)],y))
            if (np.dot(char_coordinate(y,simple_root),ty))%2 == 1:
                D.append(y)
        A.append(D)
    mys = [C[0]]
    for i in range(len(A)) :
        p = A[i]
        for j in range(len(A)) :
            q = A[j]
            if C[i] in mys:
                if len(set(p).intersection(set(q))) != 0 :
                    mys.append(C[j])
    mys = list(set(mys))
    number_of_components = len(C)//len(mys)
    if number_of_components == 1:
        print("connected")
    else:
        print("not connected : number of components are %d"%number_of_components)
    return number_of_components




def maximal_simplex_component(simple_root, fundamental_coweight, type_number):
    simple_root_sub = simple_root[1:]
    ty = Type(simple_root)[0][type_number-1]
    vtxsetp = vtx(simple_root_sub, fundamental_coweight)
    group = coset_representation(simple_root, fundamental_coweight)
    cone_label = []
    cone = []
    for i in range(len(group)) :
        ele = group[i]
        y = fundamental_coweight[0][0]
        for m in ele :
            y = tuple(reflection(simple_root[int(m)],y))
        if (np.dot(char_coordinate(y,simple_root),ty))%2 == 1 :
            cone_label.append(i)
            cone.append({i})
            break
    if len(cone_label)!=0:
            g = group[cone_label[0]]
    for i in range(len(group)) :
        ele = group[i]
        y = fundamental_coweight[0][0]
        for m in range(len(ele)) :
            y = tuple(reflection(simple_root[int(ele[m])],y))
        if (np.dot(char_coordinate(y,simple_root),ty))%2 == 1 :
            count = 0
            for x in vtxsetp :
                z = tuple(list(x))
                for m in range(len(ele)) :
                    z = tuple(reflection(simple_root[int(ele[m])],z))
                for m in range(len(g)) :
                    x = tuple(reflection(simple_root[int(g[m])],x))
                if (np.dot(char_coordinate(x,simple_root),ty))%2 == (np.dot(char_coordinate(z,simple_root),ty))%2 :
                    count +=1
            if not count==len(vtxsetp) :
                cone_label.append(i)
                g = group[cone_label[-1]]
                cone.append({i})
            else :
                cone[-1].add(i)
    A = []
    C = []
    for i in cone:
        C += list(i)
    for i in C:
        ele = group[i]
        D = []
        for j in range(len(vtxsetp)):
            y = vtxsetp[j]
            for m in ele:
                y = tuple(reflection(simple_root[int(m)],y))
            if (np.dot(char_coordinate(y,simple_root),ty))%2 == 1:
                D.append(y)
        A.append(D)
    mys = [C[0]]
    for i in range(len(A)) :
        p = A[i]
        for j in range(len(A)) :
            q = A[j]
            if C[i] in mys:
                if len(set(p).intersection(set(q))) != 0 :
                    mys.append(C[j])
    mys = list(set(mys))
    mys.sort()
    T = set()
    for i in mys:
        with open('simplices of {} of type_{} No.{} small.pkl'.format(name, type_number,i),'rb') as h:
            lst = pickle.load(h)
        T = T.union(set(lst))
    S = set()
    for x in T:
        S = S.union(set(x))
    U = set()
    for i in C:
        with open('simplices of {} of type_{} No.{} small.pkl'.format(name, type_number,i),'rb') as h:
            lst = pickle.load(h)
        for x in lst:
            if len(S.intersection(set(x))) != 0:
                U.add(x)
    T.union(U)
    not_cone = []
    for i in range(len(group)):
        if i not in C:
            not_cone.append(i)
    LST = set()
    for i in not_cone:
        with open('simplices of {} of type_{} No.{} small.pkl'.format(name, type_number,i),'rb') as h:
            lst = pickle.load(h)
        LST = LST.union(set(lst))
    V = set()
    S = set()
    for x in list(T):
        S = S.union(set(x))        
    for x in LST:
        if len(S.intersection(set(x))) != 0:
            V.add(x)
    lst_component = T.union(V)
    with open('simplices of {} of type_{} component small.pkl'.format(name, type_number), 'wb') as h:
        pickle.dump(list(lst_component), h, pickle.HIGHEST_PROTOCOL)
    print('complete.', 'list size :', len(lst_component))




# Vertex Elimination

def maximal_simplex_VE(simple_root, fundamental_coweight):
    typ = Type(simple_root)
    vtxset = vtxpre(simple_root, fundamental_coweight)
    L = vtxpre(simple_root, fundamental_coweight)
    L.sort(key = len)
    order = []
    for i in L:
        order.append(vtxset.index(i))
    type_number = 1
    for ty in typ[0]:
        start_sub = time.time()
        if number_of_components(simple_root, fundamental_coweight, type_number) == 1:
            with open('simplices of {} of type_{} small.pkl'.format(name, type_number), 'rb') as h:
                sim = pickle.load(h)
        else:
            with open('simplices of {} of type_{} component small.pkl'.format(name, type_number), 'rb') as h:
                sim = pickle.load(h)
        for i in order:
            mtx = []
            for x in vtxset[i]:
                if np.dot(ty, char_coordinate(x,simple_root))%2 == 1:
                    T = set()
                    for j in sim:
                        if x in j:
                            y = list(j)
                            y.remove(x)
                            T.add(tuple(y))
                    mtx.append(list(T))
                else:
                    mtx.append([])
            count = 0
            T = []
            for k in range(len(mtx)):
                lst = mtx[k]
                S = sc.SimplicialComplex(lst)
                c = 0
                for j in range(dim):
                    if S.betti(j) == 0:
                        c+=1
                if c == dim :
                    count += 1
                else:
                    T.append(k)
            if count != 0:
                if len(T) == 0:
                    S = set(vtxset[i])
                else:
                    S = set()
                    for k in range(len(vtxset[i])):
                        if k not in T:
                            S.add(vtxset[i][k])
                pim = set()
                for x in sim:
                    pim.add(tuple(set(x).difference(S)))
                sim = list(pim)
        if number_of_components(simple_root, fundamental_coweight, type_number) == 1:
            with open('simplices of {} of type_{} vtx elimination ver.pkl'.format(name, type_number), 'wb') as h:
                pickle.dump(sim, h, pickle.HIGHEST_PROTOCOL)
        else:
            with open('simplices of {} of type_{} component vtx elimination ver.pkl'.format(name, type_number), 'wb') as h:
                pickle.dump(sim, h, pickle.HIGHEST_PROTOCOL)
        print('type_{} complete'.format(type_number), 'list size :', len(sim))
        type_number += 1






