void thomas(const int X, double x[restrict X],
            const double a[restrict X], const double b[restrict X],
            const double c[restrict X], double scratch[restrict X]) {
    /*
     solves Ax = d, where A is a tridiagonal matrix consisting of vectors a, b, c
     X = number of equations
     x[] = initially contains the input v, and returns x. indexed from [0, ..., X - 1]
     a[] = subdiagonal, indexed from [1, ..., X - 1]
     b[] = main diagonal, indexed from [0, ..., X - 1]
     c[] = superdiagonal, indexed from [0, ..., X - 2]
     scratch[] = scratch space of length X, provided by caller, allowing a, b, c to be const
     not performed in this example: manual expensive common subexpression elimination
     */
    scratch[0] = c[0] / b[0];
    x[0] = x[0] / b[0];

    /* loop from 1 to X - 1 inclusive */
    for (int ix = 1; ix < X; ix++) {
        if (ix < X-1){
        scratch[ix] = c[ix] / (b[ix] - a[ix] * scratch[ix - 1]);
        }
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) / (b[ix] - a[ix] * scratch[ix - 1]);
    }

    /* loop from X - 2 to 0 inclusive */
    for (int ix = X - 2; ix >= 0; ix--)
        x[ix] -= scratch[ix] * x[ix + 1];
}