import numpy as np


# atom struct
atom_dtype = np.dtype([
        ('r', 'f8', (3,)), # position
        ('m', 'f8'), # mass
        ('elem', 'U2') # element
])

# bond struct
bond_dtype = np.dtype([
        ('a1', int), # atom 1 index
        ('a2', int), # atom 2 index
        ('R', int, (2,)), # unit cell index
        ('l', 'f8'), # length of bond
        ('d', 'f8', (3,)), # direction of bond
        ('k', 'f8') # strength of bond
])

def RotationMatrix(theta):
        return np.array([[np.cos(theta), np.sin(theta), 0],
                        [-np.sin(theta), np.cos(theta), 0],
                        [0, 0, 1]])


def GenerateMomentumLoop(q_points, segment_n):
        q_loop = []
        for i in range(len(q_points)):
                start = np.array(q_points[i])
                end = np.array(q_points[(i+1) % len(q_points)])

                for j in range(segment_n):
                        q_loop.append(start + j * (end-start)/segment_n)
        return q_loop


def CreateBonds(unit_cell, neighboring_unit_cells, lattice_x, lattice_y, bonding_params):
        d_cutoff = bonding_params[0]
        inter_layer_strength = bonding_params[1]
        inter_layer_decay = bonding_params[2]
        intra_layer_strength = bonding_params[3]
        intra_layer_decay = bonding_params[4]

        bonds = []
        for a1 in range(len(unit_cell)):
                for a2 in range(len(unit_cell)):
                        for R in neighboring_unit_cells:
                                vec = unit_cell[a2]['r'] + lattice_x * R[0] + lattice_y * R[1] - unit_cell[a1]['r']
                                if np.linalg.norm(vec) <= d_cutoff and np.linalg.norm(vec) != 0:
                                        # Inter layer
                                        if vec[2] > .1:
                                                bonds.append(
                                                        (a1, # Atom 1 index
                                                         a2, # Atom 2 index
                                                         R, # Unit Cell index
                                                         np.linalg.norm(vec), # length of bond
                                                         vec/np.linalg.norm(vec), # Direction of bond
                                                         inter_layer_strength * np.exp(-inter_layer_decay * np.linalg.norm(vec)) # Strength of bond
                                                        )
                                                )
                                        # Intra layer
                                        else:
                                                bonds.append(
                                                        (a1, # Atom 1 index
                                                         a2, # Atom 2 index
                                                         R, # Unit Cell index
                                                         np.linalg.norm(vec), # length of bond
                                                         vec/np.linalg.norm(vec), # Direction of bond
                                                         intra_layer_strength * np.exp(-intra_layer_decay * np.linalg.norm(vec))) # Strength of bond
                                                )

        return bonds


def DynamicsMatrix(bonds, unit_cell, lattice_x, lattice_y, q):
        dynamics_matrix = np.zeros((3 * len(unit_cell), 3* len(unit_cell)), dtype = complex)

        for b_i, bond in enumerate(bonds):
                a1 = bond[0]
                a2 = bond[1]
                R = bond[2]
                l = bond[3]
                vec = bond[4]
                d3x3 = bond[5] * np.outer(vec,vec)

                # Only grabbing z component
                #d3x3 = bond[5]  * np.array([[0,0,0],[0,0,0],[0,0,vec[2]*vec[2]]])
                m = unit_cell[a1]['m']

                dynamics_matrix[3*a1:3*a1+3, 3*a2:3*a2+3] -= d3x3 * np.exp(1j * np.dot(q, R[0] * lattice_x + R[1] * lattice_y))/m # Off diagonal
                dynamics_matrix[3*a1:3*a1+3, 3*a1:3*a1+3] += d3x3/m # Diagonal
        return dynamics_matrix
