import numpy as np
import Base


### ONLY WORKS FOR SQUARE LATTICE
def GenerateMoireParams(legx, legy, h):
        """
        Create a Moire pattern based on the given legx and legy values.
        Parameters:
        legx (int): x leg of Moire pattern.
        legy (int): y leg of Moire pattern.
        Returns:
        tuple: Contains twist angle, bottom, top, and Moire lattice vectors.
        """

        # Check if the Moire pattern is commensurate
        if (int(np.sqrt(legx**2 + legy**2)) != np.sqrt(legx**2 + legy**2)):
                print("Warning: Moire pattern is not commensurate")


        # Twist Angle   
        theta = np.atan(legx/legy)
        R = Base.RotationMatrix(theta)

        # Bottom lattice vectors
        ebx = np.array([1, 0, 0]);
        eby = np.array([0, 1, 0]);


        # Top lattice vectors
        etx = np.dot(R, ebx)
        ety = np.dot(R, eby)

        # Moire Lattice Vectors 
        eMx = np.sqrt(legx**2 + legy**2) * etx;
        eMy = np.sqrt(legx**2 + legy**2) * ety;

        # Moire inverse lattice vectors
        GMx = 2 * np.pi * Base.RotationMatrix(np.pi/2) @ eMy / (eMx @ Base.RotationMatrix(np.pi/2) @ eMy)
        GMy = 2 * np.pi * Base.RotationMatrix(np.pi/2) @ eMx / (eMy @ Base.RotationMatrix(np.pi/2) @ eMx)


        return theta, h, ebx, eby, etx, ety, eMx, eMy, GMx, GMy



def GenerateMoireUnitCell(unit_cell, moire_params, delta = .01):
        """
        unit_cell : specified as list of [x,y,z,m] for each atom
        moire_params : twist, height, lattice vectors (top, bottom, moire)
        """
        # Unpack params
        twist, h, ebx, eby, etx, ety, eMx, eMy = moire_params

        bottom_unit_cell = unit_cell
        top_unit_cell = np.array([(Base.RotationMatrix(twist) @ a['r'], a['m'], a['elem']) for a in unit_cell], dtype=Base.atom_dtype)

        n = int(eMx@eMx + 1)

        moire_unit_cell = []

        for i in range(-n, n, 1):
                for j in range(-n, n, 1):
                        # Bottom layer
                        r = i * ebx + j * eby
                        # If within moire lattice vectoprs
                        if r @ eMx / (eMx @ eMx) <= 1-delta and r @ eMx / (eMx @ eMx) > -delta and r @ eMy / (eMy @ eMy) <= 1-delta and r @ eMy / (eMy @ eMy) > -delta:
                                # For atom
                                for a in bottom_unit_cell:
                                        moire_unit_cell.append((a['r'] + r, a['m'], a['elem']))
                        # Top layer
                        r = i * etx + j * ety
                        # If within moire
                        if r @ eMx / (eMx @ eMx) <= 1-delta and r @ eMx / (eMx @ eMx) > -delta and r @ eMy / (eMy @ eMy) <= 1-delta and r @ eMy / (eMy @ eMy) > -delta:
                                # For atom in top unit cell
                                for a in top_unit_cell:
                                        moire_unit_cell.append((a['r'] + r + np.array([0,0,h]), a['m'], a['elem']))

        return np.array(moire_unit_cell, dtype = Base.atom_dtype)
