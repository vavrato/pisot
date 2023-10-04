proof.number_field(False) #GRH is assumed, makes computation of fundamental units run much faster
global verbose     #set 1 for text output/debug
global precision
precision = 250    #needs to be very big sometimes, for example for x^4-101*x^3+5*x^2-101*x+1 the precision must be around 250; otherwise sage makes some errors

def msg(text):     #text output for verbose mode
    if verbose:
        print(text)

def inicialize(f,x):
    ''' Initialization of the fields, matrices, ...

    Returns the Log-matrix, the embedding corresponding to the real generator closest to x, the corresponding number field, its fundamental system of units

    '''
    NF.<g> = NumberField(f,embedding = x)
    print(g)
    msg('Computing the units...')
    units = NF.units()
    msg('DONE!')
    msg('-------------------')
    emb = g.complex_embeddings(precision)
    ran = f.degree()-1
    indices = [(j,emb[j].is_real()) for j in range(ran) if emb[j].conjugate() != emb[j+1]] + [(ran,emb[ran].is_real())]    #indices of embeddings when each first of a complex pair is omitted; together with an information if the embedding is real

    emb_nopairs= [emb[j] for j,r in indices]
    k = [abs(i-g.n(precision))<10^(-30) for i in emb_nopairs].index(1) # k is the embedding corresponding to the generator, BEWARE: the 'e-30' precision in the term 10^(-30) must be somewhat smaller than the global precision

    cols = []
    for u in units:
        cols.append([(2^(1-r))*log(abs(u.complex_embeddings(precision)[j])) for j,r in indices])
    A = column_matrix(cols)

    msg('The roots of the polynomial are (not showing the complex conjugates): ')
    if verbose:
        for i in range(len(emb_nopairs)):
            print(str(i) + ': ' + str(emb_nopairs[i]))

    msg('\n' + 'Your field is determined by the ' +str(k)+'-th root in the previous list (indexing from zero)')
    msg('-------------------')

    msg('Units of the field are: ' + str(units))
    msg('-------------------')

    msg('The Log matrix is: ')
    msg(str(A))
    msg('-------------------')

    return [A,k,NF,units]


def FINDMIN(A,k,N,units):
    ''' Finds minimal U-number in the field N specified by the k-th embedding (see inicialize())

    returns the minimal U-number of the field, and its vector of exponents

    '''
    deg = N.absolute_degree()
    if (deg == 2) or (deg % 2 == 1):
        delta = log(1.32)
    else:
        delta = (1/4)*(log(log(deg))/log(deg))^3
    p = MixedIntegerLinearProgram(maximization = False)    #beware: glpk solver, although being quite fast, sometimes returns a non-admissible solution

    B = copy(A)

    B.set_row_to_multiple_of_row(k,k,-1)

    x = p.new_variable(nonnegative=False, integer = True)
    obj = [-B[k][i]*x[i] for i in range(len(B[0]))]    #objective function set to be k-th row of original matrix times x
    p.set_objective(sum(obj))
    c = [0] * B.nrows()
    c[k] = -delta
#    c[k] = - 0.97242365020

    p.add_constraint(B * x <= c)
    if verbose ==1:
        print('FINDMIN solves the following: ' + '\n')
        p.show()
        print('-------------------')
    p.solve()
    e = vector(p.get_values(x).values())

    unumber = prod([units[i]^e[i] for i in range(len(e))])
    unumber = unumber if unumber.n() > 0 else -unumber

    msg('The optimal value of FINDMIN is achieved at: ' + str(e))
    msg('The U-number output of FINDMIN is: ' + str(unumber.n()) + ' with minimal polynomial ' + str(unumber.minpoly()))

    return (unumber,e) if unumber.n() > 0 else (-unumber,e)

def CUTEDGE(A,k,N,units):
    ''' Implementation of CUTEDGE

    returns the Pisot generator of the number field N of minimal height

    '''
    deg = N.absolute_degree()
    if (deg == 2) or (deg % 2 == 1):
        delta = log(1.32)
    else:
        delta = (1/4)*(log(log(deg))/log(deg))^3
    M = copy(A)
    out = FINDMIN(M,k,N,units)
    unumber = out[0]
    e = out[1]

#    unumber = prod([units[i]^e[i] for i in xrange(len(e))])
#    unumber = unumber if unumber.N() > 0 else -unumber

#    msg('The optimal value of FINDMIN is achieved at: ' + str(e))
#    msg('The U-number output of FINDMIN is: ' + str(unumber.N()) + ' with minimal polynomial ' + str(unumber.minpoly()))


    if 1.0 not in map(abs,unumber.complex_embeddings(prec=300)):
        msg('\t'+'* and it is (complex) Pisot with conjugates:')
        msg('\t'+str(unumber.complex_embeddings()))
        if unumber in RR:
            msg('!!!WARNING!!!')
        return unumber if unumber.n() > 0 else -unumber
    else:
        msg('\t'+'* and it is Salem with conjugates:')
        msg('\t'+str(unumber.complex_embeddings()))
        msg('\t'+'* and absolute values of the conjugates:')
        msg('\t'+str(map(abs,unumber.complex_embeddings())))

        j = 0     # j will be the the non-identic real embedding
        for i in range(A.nrows()):
            x = A[i]*e
            if (x<0) & (abs(x)>10^(-50)):
                j = i
                msg('\t'+'* and the other real embedding is the ' + str(j) + '-th one. (In the first list).')
                break

    msg('-------------------')

    B = copy(A)
    B.add_multiple_of_row(k,j,1)
    B.set_row_to_multiple_of_row(k,k,-1)    # now we have C^k = -A^k - A^j
    c = [0] * B.nrows()
    c[k] = -delta


    q = MixedIntegerLinearProgram(maximization = False)
    x = q.new_variable(nonnegative=False, integer = True)
    obj = [A[k][i]*x[i] for i in range(len(A[0]))]
    q.set_objective(sum(obj))
    q.add_constraint(B * x <= c)
    if verbose == 1:
        print('The second run of FINDMIN solves the following: ' + '\n')
        q.show()
        print('-------------------')

    q.solve()
    e = vector(q.get_values(x).values())

    pisot = prod([units[i]^e[i] for i in range(len(e))])
    pisot = pisot if pisot.n() > 0 else -pisot

    msg('The optimal value of CUTEDGE is achieved at: ' + str(e))
    msg('The smallest Pisot unit in your field is approximately ' + str(pisot.n()) + ' with minimal polynomial ' + str(pisot.minpoly()))

    return pisot
###################
verbose = 1
P.<x> = PolynomialRing(ZZ)
f = x^2-5*x+1 # defining polynomial for the field
emb = 1 # we choose an embedding
i = inicialize(f,emb)
CUTEDGE(i[0],i[1],i[2],i[3]).n(precision)