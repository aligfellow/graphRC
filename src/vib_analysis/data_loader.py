import json
from importlib import resources

BOHR_TO_ANGSTROM = 0.52917721054 
 
def load_vdw_radii():
    """Load VDW radii data from package data directory"""
    data_path = resources.files('vib_analysis.data')
    vdw_file = data_path / 'vdw_radii.json'
    with vdw_file.open('r') as f:
        data = json.load(f)
    
    return {element: data[element]['vdw_radius'] * BOHR_TO_ANGSTROM
            for element in data if not element.startswith('_')}

def load_expected_valences():
    """Load expected valences from JSON file"""
    try:
        data_path = resources.files('vib_analysis.data')
        valence_file = data_path / 'expected_valences.json'
        
        with open(valence_file, 'r') as f:
            data = json.load(f)
        
        return {element: data[element] for element in data if not element.startswith('_')}

    except Exception as e:
        print(f"Error loading valences file: {e}, using minimal fallback")
        return {
            'H': [1], 'C': [4], 'N': [3, 5], 'O': [2], 'F': [1],
            'P': [3,5], 'S': [2,4,6], 'Cl': [1], 'Br': [1], 'I': [1,2,3],
            'Li': [1], 'Na': [1], 'K': [1], 'Mg': [2], 'Ca': [2],
            'Zn': [2], 'Fe': [2, 3], 'Cu': [1, 2], 'Mn': [2, 3, 4], 'Ru': [2, 3, 4]
        }
