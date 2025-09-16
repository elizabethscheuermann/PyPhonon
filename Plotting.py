import matplotlib.pyplot as plt
import numpy as np

def PlotUnitCell(unit_cell):
        fig = plt.figure()

        ax = fig.add_subplot(111,projection='3d')

        for atom in unit_cell:
                ax.scatter(atom['r'][0], atom['r'][1], atom['r'][2], color = atom_colors[atom['elem']])


        plt.show()

def PlotBonds(bonds, unit_cell, lattice_x, lattice_y):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        for bond in bonds:
                x1 = unit_cell[bond[0]]['r'][0] # X Position of atom 1
                y1 = unit_cell[bond[0]]['r'][1] # Y position of atom 1
                z1 = unit_cell[bond[0]]['r'][2] # Z position of atom 1
                x2 = unit_cell[bond[1]]['r'][0] # X Position of atom 2
                y2 = unit_cell[bond[1]]['r'][1] # Y Position of atom 2
                z2 = unit_cell[bond[1]]['r'][2] # Z position of atom 2
                R = bond[2][0] * lattice_x + bond[2][1] * lattice_y

                ax.plot3D([x1, x2 + R[0]], [y1, y2 + R[1]], [z1, z2 + R[2]], color = 'black')

        plt.show()



def PlotBandStructure(bands):

        fig, ax = plt.subplots()
        bands = np.real(np.array(bands))

        for i in range(len(bands[0])):
                ax.plot(bands[:, i], linestyle = '', marker='.', color = 'black')
        plt.show()
