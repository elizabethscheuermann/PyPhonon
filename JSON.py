import json


def LoadParamsFromJSON(filepath):
	with open(filepath, 'r') as file:
		data = json.load(file)

	return data


data = LoadParamsFromJSON("test.json")


print(data['lattice_x'])


