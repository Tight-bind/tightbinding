from libc.math cimport sqrt

cdef inline norm(double [:] tdvect):
    return (sqrt(tdvect[0]*tdvect[0] + tdvect[1]*tdvect[1] +
            tdvect[2]*tdvect[2]))
