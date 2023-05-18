import numpy as np


def normalize(p):
    return np.array(p) /np.sqrt(np.dot(p, p))


p1 = normalize([1, 2, 3])
p2 = normalize([4, 5, 6])
p3 = normalize([7, 8, 1])


p1 = [0.15329, 0.320602, -0.934728]
p2 = [0.734215, -0.249426, 0.631439]
p3 = [-0.876535, -0.431829, -0.212626]


rho = np.sqrt(2 * (1+np.dot(p1,p2)) * (1+np.dot(p2,p3)) * (1+np.dot(p3,p1)))
print(np.dot(p1,p2), np.dot(p2,p3), np.dot(p3,p1))
print(rho)
a = (1 + np.dot(p1, p2) + np.dot(p2, p3) + np.dot(p3, p1)) / rho
b = (np.dot(p1, np.cross(p3, p2))) / rho

s, c = np.arccos(a), np.arcsin(b)
print(s, " // ", a)
print(c, " // ", b)



# for n in range(10):
#     p1 = normalize(np.random.rand(3))
#     p2 = normalize(np.random.rand(3))
#     p3 = normalize(np.random.rand(3))
#
#     rho = np.sqrt(2 * (1+np.dot(p1,p2)) * (1+np.dot(p2,p3)) * (1+np.dot(p3,p1)))
#     a = (1 + np.dot(p1, p2) + np.dot(p2, p3) + np.dot(p3, p1)) / rho
#     b = (np.dot(p1, np.cross(p3, p2))) / rho
#
#     s, c = np.arccos(a), np.arcsin(b)
#     if not np.isclose(abs(s), abs(c), atol=1e-32):
#         print(s, " // ", c)
