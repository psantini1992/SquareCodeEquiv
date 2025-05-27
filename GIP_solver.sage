from collections import defaultdict

#Computes the graph associated with generator matrix
#Works only if G has trivial hull, i.e., rank(G*G^T) = k
def graph_from_generator_matrix(G):
    X = G*G.transpose()
    return G.transpose()*X^-1*G

###########################################################

def count_in_Fq(n, Fq_list, a):
    vals = vector(ZZ, len(Fq_list))
    for i in range(n):
        j = Fq_list.index(a[i])
        vals[j] += 1
    return vals

#Return decision (either True or False) and candidates for images of indices
def solve_GIP(Fq, A_1, A_2):
    
    Fq_list = Fq.list()
    n = A_1.ncols()
    
    #Now, check the multiset formed by all elements: they must be the same
#    multiset_1 = []; multiset_2 = []
#    for i in range(0,n-1):
#        for j in range(i+1,n):
#            multiset_1.append(A_1[i,j])
#            multiset_2.append(A_2[i,j])
    
#    if sorted(multiset_1) != sorted(multiset_2):
#        return False, []
    
    
    #Build a list with all entries of a row, hashed with the element in the diagonal
    #We're gonna use a Python dictionary because it makes finding collisions easier
    #We use the hash result as the key for the dictionary
    rows_labels_1 = defaultdict(list)
    for i in range(n):
        vals = count_in_Fq(n, Fq_list, [A_1[i,j] for j in range(n) ]) #all elements in the row
        label = str(vals) #hash to simplify storing the element
        label = str(vals)+str(A_1[i,i]) #concatenate with element in main diagonal and hash again
        if rows_labels_1.get(label) == None:
            rows_labels_1[label] = [i]
        else:
            rows_labels_1[label].append(i)

    #Repeat for second graph
    rows_labels_2 = defaultdict(list)
    for i in range(n):
        vals = count_in_Fq(n, Fq_list, [A_2[i,j] for j in range(n) ])
        label = str(vals)
        label = str(vals)+str(A_2[i,i])
        if rows_labels_2.get(label) == None:
            rows_labels_2[label] = [i]
        else:
            rows_labels_2[label].append(i)
    
    
    #build collisions, check also for isomorphism not existing
    rows_colls = []
    for label in rows_labels_1.keys():
        indices_1 = rows_labels_1.get(label)
        indices_2 = rows_labels_2.get(label)
        
        #If number of collisions is not the same, return False
        if indices_1 == None:
            len_indices_1 = 0
        else:
            len_indices_1 = len(indices_1)

        if indices_2 == None:
            len_indices_2 = 0
        else:
            len_indices_2 = len(indices_2)

        if len_indices_1 != len_indices_2:
            return False, []
        rows_colls.append([indices_1, indices_2])
    
    #All test have been successful, return ok and found collisions
    return True, rows_colls
    
####################################################

#This function assumes every row from A_1 has a unique collision from A_2
def simple_solve_GIP(Fq, A_1, A_2):

    Fq_list = Fq.list()
    n = A_1.ncols()


    #Build a list with all entries of a row, hashed with the element in the diagonal
    #We're gonna use a Python dictionary because it makes finding collisions easier
    #We use the hash result as the key for the dictionary
    rows_labels_1 = defaultdict(list)
    for i in range(n):
        a_i = A_1[i].list()
        a_i.sort()
        rows_labels_1[str(a_i)] = i

    #Compute multisets for second graph, check for collisions
    rows_colls = []
    for j in range(n):
        b_j = A_2[j].list()
        b_j.sort()

        #find index of collision
        i = rows_labels_1.get(str(b_j))
        rows_colls.append([i,j])


    #All test have been successful, return found collisions
    return rows_colls

####################################################
