import math
import pylinear.matrices as num
import pylinear.linear_algebra as la

n1 = num.array([-0.48541,-0.352671])
n2 = num.array([-0.525784,-0.289052])
n3 = num.array([-0.594783,-0.377461])

pre_g = num.array([[1,1,1],[n1[0], n2[0], n3[0]], [n1[1], n2[1], n3[1]]])
pre_g_inv = la.inverse(pre_g)
t_matrix = num.array([[n2[0]-n1[0], n3[0]-n1[0]],[n2[1]-n1[1],n3[1]-n1[1]]])
g_part_2=num.array([[0,0],[1,0],[0,1]])
g = num.matrixmultiply(pre_g_inv, g_part_2)
m = 0.5 * num.matrixmultiply(g, num.transpose(g)) * math.fabs(la.determinant(t_matrix))
print m
