def square_code(G):
    '''
    Takes as input generator G for linear codes, returns generator for square code in systematic form
    '''
    k = G.nrows(); n = G.ncols(); Fq = G[0,0].parent()
    new_G = matrix(Fq,k+k*(k-1)/2,n)
    row_index = 0
    for i in range(k):
        for j in range(i,k):
            for ell in range(n):
                new_G[row_index,ell] = G[i,ell]*G[j,ell]
            row_index+=1

    #NOTE: the output matrix may have less than k(k+1)/2 rows
    return codes.LinearCode(new_G).generator_matrix()


def hull(Fq, G):
    '''
    Compute the hull of the code generated by input G, defined over the input finite field Fq
    '''
    H = codes.LinearCode(G).parity_check_matrix()
    r = H.nrows()

    #To compute hull, stack G on top of H and compute right kernel
    A = block_matrix(Fq,2,1,[G,H])
    B = A.right_kernel()

    return matrix(B.basis()) #return a matrix

def hull_dimension(G):
	'''
	Computes dimension of hull of the code generated by the input matrix
	'''
	k = G.nrows()
	hull_dim = k - rank(G*G.transpose())

	return hull_dim
