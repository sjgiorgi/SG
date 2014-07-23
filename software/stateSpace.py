#!/usr/bin/env python

class stateSpace(object):
    def __init__(self, *args):
        if len(args) == 5:
            (A, B, C, D, Ts) = args
            if Ts == 0:
                Ts = None
        elif len(args) == 4:
            (A, B, C, D) == args
            Ts = None
        elif len(args) == 1:
            print type(args[0])
            if not isinstance(args[0], stateSpace):
                raise TypeError("Argument must be stateSpace. Received %s." % type(args[0]))
            A = args[0].A
            B = args[0].B
            C = args[0].C
            D = args[0].D
            try:
                Ts = args[0].Ts
                if Ts == 0:
                    Ts = None
            except NameError:
                Ts = None
        else:
            print "Must be either 1, 4, or 5 arguements"
        # check / convert arguements to matrices
        inputs = [A,B,C,D]
        for i in range(len(inputs)):
            if not isinstance(inputs[i],sage.matrix.matrix_generic_dense.Matrix_generic_dense):
                try:
                    inputs[i] = matrix(inputs[i])
                except TypeError:("Arguments must be of the form: sage.matrix.matrix_generic_dense.Matrix_generic_dense")
        [A,B,C,D] = inputs

        # check matrix dimensions
        if A.nrows() != A.ncols():
            raise TypeError("A must be square.")
        if B.nrows() != A.nrows():
            raise TypeError("Number of rows of B must match A.")
        if C.ncols() != A.nrows():
            raise TypeError("Number of cols of C must match A.")
        if D.nrows() != C.nrows():
            raise TypeError("Number of rows of D must match C")
        if D.ncols() != B.ncols():
            raise TypeError("Number of cols of D must match B")
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.Ts = Ts

    # series connection
    def seriesConnection(self, secondSS):
        # check if secondSS is a stateSpace object
        if not isinstance(secondSS, stateSpace) and not isinstance(secondSS, transFunct):
            raise TypeError("Argument must be stateSpace or transFunct. Received %s." % type(secondSS))
        # convert to SS if necessary
        if isinstance(secondSS, transFunct):
            secondSS = tf2ss(secondSS)
        return self

        # convert if necessary

    # parallel connection
    def parallelConnection(self, secondSS):
        # check if secondSS is a stateSpace object
        if not isinstance(secondSS, stateSpace) and not isinstance(secondSS, transFunct):
            TypeError("Argument must be stateSpace or transFunct. Received %s." % type(secondSS))
        # convert to SS if necessary
        if isinstance(secondSS, transFunct):
            secondSS = tf2ss(secondSS)
        return self

    def isControllable(self):
        controlMatrix = self.controllabilityMatrix()
        return self.A.rank() == controlMatrix.rank()

    def controllabilityMatrix(self):
        controlMatrix = self.B
        for i in range(1,self.A.ncols()):
            controlMatrix = block_matrix([[ controlMatrix, (self.A^i)*self.B ]])
        return controlMatrix

    def numOfUncontrolableState(self):
        controlMatrix = self.controllabilityMatrix()
        return self.A.nrows()-controlMatrix.rank()

    def isObservable(self):
        observeMatrix = self.observabilityMatrix()
        return self.A.rank() == observeMatrix.rank()

    def observabilityMatrix(self):
        observeMatrix = self.C
        for i in range(1,self.A.ncols()):
            observeMatrix = block_matrix([[observeMatrix],[self.C*self.A^i ]])
        return observeMatrix

    def numOfUnobservableState(self):
        observeMatrix = self.observabilityMatrix()
        return self.A.nrows()-observeMatrix.rank()

    # evaluate state space at at a list of frequencies
    def evalFreq(self, freq):
        return self

    # returns poles of state space
    def poles(self):
        return self

    # returns zeros of state space
    def zeros(self):
        return self

    def isDiscrete(self):
        return self.Ts != None

    def isContinuous(self):
        return self.Ts == None

class transFunct(object):
    def __init__(self, *args):
        if len(args) == 3:
            (num, den, Ts) = args
        elif len(args) == 2:
            (num, den) = args
            Ts = 0
        elif len(args) == 1:
            if not isinstance(args[0], transFunct):
                raise TypeError("Argument must be transFunct. Received %s." % type(args[0]))
            num = args[0].num
            den = args[0].den
            try:
                Ts = args[0].Ts
            except NameError:
                Ts = 0;
        else:
            raise TypeError("Must be either 1, 2, or 3 arguments")

        # check / convert arguements to numpy arrays
        inputs = [num, den]
        for i in range(len(inputs)):
            if not isinstance(inputs[i],numpy.ndarray):
                try:
                    inputs[i] = numpy.array(inputs[i])
                except TypeError:("Arguments must be of the form: numpy.ndarray")
        [num, den] = inputs

        # check / convert sample time
        if not isinstance(Ts, sage.rings.real_mpfr.RealNumber):
            try:
                Ts = RR(Ts)
            except TypeError:("Arguments must be of the type: sage.rings.real_mpfr.RealNumber")
        self.num = num
        self.den = den
        self.Ts = Ts

    # evaluate transfer function at a list of frequencies
    def evalFreq(self, freq):
        return 

    # returns poles of transfer function
    def poles(self):
        return 

    # returns zeros of transfer function
    def zeros(self):
        return

    def isDiscrete(self):
        return self.Ts != None

    def isContinuous(self):
        return self.Ts == None

def ss2tf(*args):
    if not isinstance(args[0], stateSpace):
        try:
            ss = stateSpace(args[0])
        except TypeError:("Single argument must be of type stateSpace. Received %s." % type(args[0]))
    else:
        ss = args[0]
    H = ss.C*(x*matrix.identity(ss.A.nrows())-ss.A)^(-1)*ss.B+ss.D
    H = H[0,0].simplify_rational()
    return transFunct(H.numerator(),H.denominator(),ss.Ts)