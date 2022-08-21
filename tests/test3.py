import sys
sys.path.append('./src')
import numpy as np
from chemsynthcalc import ChemicalReaction
from chemsynthcalc.reaction_balance import Balancer
from chemsynthcalc.chemutils import intify_coefficients

r1  = "PtNH3BrNO3+CuNH3KNO3+BeCO3=C44H50O14.98+Cu(NO3)2+PtO3+Br1.99NO2+K1.97O+BeO+HNO3"
r2  = "PtNH3BrNO3+(Fe(CN)3)4(Fe(CN)2)3+C44H50O15+Cu(NO3)2+K2BeO2=CuNH3KNO3+BeCO3+K3.97Fe(CN)6+PtO2+Br1.96NO2+HNO3"
reaction = ChemicalReaction(r2)
balance = Balancer(reaction.reactant_matrix, reaction.product_matrix, reaction.rounding_order)
#coefs = balance.calculate_coefficients_Risteski()
#print(coefs)
#print(balance.intify_coefficients(coefs, 100000))
#print(intify_coefficients(coefs))

reaction_matrix = np.array([
     [ 1,     0,    0,    0,    0,    0,    0,    0,    1,    0,    0  ],
     [ 2,   18,    0,    2,    0,    2,    0,    6,    0,    1,    1  ],
     [ 3,    0,   50,    0,    0,    3,    0,    0,    0,    0,    1  ],
     [ 1,    0,    0,    0,    0,    0,    0,    0,    0,    1.96,  0  ],
     [ 3,    0,   15,    6,    2,    3,    3,    0,    2,    2,    3  ],
     [ 0,    7,    0,    0,    0,    0,    0,    1,    0,    0,    0  ],
     [ 0,   18,   44,    0,    0,    0,    1,    6,    0,    0,    0  ],
     [ 0,    0,    0,    1,    0,    1,    0,    0,    0,    0,    0  ],
     [ 0,    0,    0,    0,    2,    1,    0,    3.97,  0,    0,    0  ],
     [ 0,    0,    0,    0,    1,    0,    1,    0,    0,    0,    0  ]
     ])

reaction_matrix[[1, 2]] = reaction_matrix[[2, 1]]
reactant_matrix = reaction_matrix[:, :5]
product_matrix = reaction_matrix[:, 5:]
MP_inverse = np.linalg.pinv(reactant_matrix)
g_matrix = np.identity(reaction_matrix.shape[0])-reactant_matrix@MP_inverse
g_matrix = g_matrix@product_matrix
y_multiply = np.linalg.pinv(g_matrix)@g_matrix
y_vector = (np.identity(y_multiply.shape[1])-y_multiply).dot(np.ones(y_multiply.shape[1]))
x_multiply = MP_inverse@reactant_matrix
x_multiply = (np.identity(x_multiply.shape[1])-x_multiply) + MP_inverse@product_matrix@y_vector.T
x_vector = x_multiply[0].T
coefs = np.squeeze(np.asarray(np.hstack((x_vector, y_vector)))).tolist()
#print(reaction.matrix)

print(coefs)
print(intify_coefficients(coefs))

print(np.equal(np.sum(reaction_matrix), np.sum(reaction.matrix)))