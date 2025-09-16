import numpy as np
import argparse
import json
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

import Base
import Moire
import Plotting


# Atom colors
atom_colors = {'H' : 'blue',
		'C' : 'black'
}

# Arg parser
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--filepath')
parser.add_argument('-d', '--distancemax')
parser.add_argument('-in', '--inlayer')
parser.add_argument('-out', '--outoflayer')
parser.add_argument('-height', '--height')
parser.add_argument('-legx', '--legx')
parser.add_argument('-legy', '--legy')


args = vars(parser.parse_args())

def LoadParamsFromJSON(filepath):
	with open(filepath, 'r') as file:
		data = json.load(file)
	
	return data


def main():	
	# Get Command Line args
	filepath = args['filepath']
	
	# Load Params from JSON
	params = LoadParamsFromJSON(filepath)

	d_cutoff = params['bond_distance_cutoff']
	h = params['layer_spacing']
	intra_strength = params['in_layer_coupling_strength']
	inter_strength = params['out_of_layer_coupling_strength']
	intra_decay = params['in_layer_coupling_decay']
	inter_decay = params['out_of_layer_coupling_decay']
	legx = params['moire_leg_x']
	legy = params['moire_leg_y']
	q_n = params['n_momentum_steps']
	momentum_points = params['momentum_loop']	

	
	# Get Moire Params
	twist, height, ebx, eby, etx, ety, eMx, eMy, GMx, GMy = Moire.GenerateMoireParams(legx, legy, h)
	
	# Create Moire Unit Cell
	moire_unit_cell = Moire.GenerateMoireUnitCell(np.array([
						((0.0, 0.0, 0.0), 1.0,'H')],
						 dtype = Base.atom_dtype),
						(twist, height, ebx, eby, etx, ety, eMx, eMy))

	# Create Momentum Loop
	points = []
	for point in momentum_points:
		if point == "Gamma":
			points.append(np.array([0,0,0]))
		elif point == "X":
			points.append(np.pi * eMx / (eMx @ eMx))
		elif point == "M":
			points.append(np.pi * eMx / (eMx @ eMx) + np.pi * eMy / (eMy @ eMy))

	q_loop = Base.GenerateMomentumLoop(points, q_n)	
	
	# Set Neighboring unit cells
	neighboring_unit_cells = [
		[-2,-2, 0],
		[-2, -1,0],
		[-2,0,0],
		[-2,1,0],
		[-2,2,0],
		[-1, -2, 0],
		[-1, -1, 0],
		[-1, 0, 0],
		[-1, 1, 0],
		[-1, 2, 0],
		[0, -2, 0],
		[0, -1, 0],
		[0, 0, 0],
		[0, 1, 0],
		[0, 2, 0],
		[1, -2, 0],
		[1, -1, 0],
		[1, 0, 0],
		[1, 1, 0],
		[1, 2, 0],
		[2, -2, 0],
		[2, -1, 0],
		[2, 0, 0],
		[2, 1, 0],
		[2, 2, 0]
	]	

	# Create Plot
	fig = plt.figure(figsize = (2,10))
	ax = fig.add_subplot(111)
	ax.set_xticks([0, len(q_loop)/3, 2*len(q_loop)/3, len(q_loop)])
	ax.set_xticklabels([r"$\Gamma$", r"$X$", r"$M$", r"$\Gamma$"])
	# Create Bonds
	bonds = Base.CreateBonds(moire_unit_cell, neighboring_unit_cells, eMx, eMy, [d_cutoff, inter_strength, inter_decay, intra_strength, intra_decay])
	ax.set_title(r"$A_{in}$ = " + str(intra_strength)
			+ r"$, \lambda_{in} = $" + str(intra_decay)
			+ r", $A_{out} = $" + str(inter_strength) 
			+ r", $\lambda_{out} = $" + str(inter_decay) 
			+ r", $h = $" + str(h) 
			+ r", $d_{max} = $" + str(d_cutoff))	

	bands = []
	
	# Momentum Loop
	for q_i, q in enumerate(q_loop):
		os.system('clear')
		print("Simulation Parameters:")
		print("Intra Strength: ", intra_strength)
		print("Intra Decay: ", intra_decay)
		print("Inter Strength: ", inter_strength)
		print("Inter Decay: ", inter_decay)
		print("Layer Spacing: ", h)
		print("Bond Cutoff Distance: ", d_cutoff)
		print("Twist Angle (degrees): ", twist*180/np.pi)
		print("Moire Leg X: ", legx)
		print("Moire Leg Y: ", legy)


		print(100*(q_i+1)/len(q_loop), "%")

		print("Creating Dynamics Matrix...")
		# Dynamics Matrix
		dm = Base.DynamicsMatrix(bonds, moire_unit_cell, eMx, eMy, q)
		print("Finding Eigenvalues...")
		# Eigvals
		vals, vecs = np.linalg.eig(dm)
		bands.append(np.sqrt(vals))
			
	# Plotting
	print("Plotting...")
	bands = np.real(np.array(bands))
	for i in range(len(bands[0])):
		ax.plot(bands[:,i], linestyle = '', marker = '.', color = "blue")

	
	plt.show()
				
if __name__ == "__main__":
	main()
