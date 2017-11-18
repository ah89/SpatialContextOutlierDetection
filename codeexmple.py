from cvxopt import matrix, solvers

Q = matrix(0.0,(3,3))
Q[::4] = 2.0
r = matrix(0.0, (3,1))

A = matrix([[1,0],[1,1],[1,0]],(2,3),tc='d')
b = matrix([6.0, 0.0])
C = matrix([1.0, 0.0, 0.0],(1,3))
d = matrix(1.0)



sol=solvers.qp(Q, r, C, d, A, b)
print("\nSolution =\n" + str(sol['x']))