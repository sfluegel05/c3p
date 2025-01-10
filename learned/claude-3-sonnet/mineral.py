"""
Classifies: CHEBI:46662 mineral
"""
"""
Classifies: CHEBI:33686 mineral
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mineral(smiles: str):
    """
    Determines if a molecule is a mineral based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a mineral, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all atoms
    atoms = mol.GetAtoms()
    if len(atoms) == 0:
        return False, "No atoms found"
    
    # Count different types of atoms
    metal_ions = set()
    has_metal = False
    organic_atoms = 0
    inorganic_atoms = 0
    charged_atoms = 0
    
    # Define common metal atoms in minerals
    metals = {11: 'Na', 12: 'Mg', 13: 'Al', 19: 'K', 20: 'Ca', 24: 'Cr',
              25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn',
              38: 'Sr', 47: 'Ag', 48: 'Cd', 56: 'Ba', 57: 'La', 78: 'Pt',
              79: 'Au', 80: 'Hg', 82: 'Pb', 83: 'Bi', 55: 'Cs', 51: 'Sb',
              46: 'Pd', 33: 'As'}
    
    # Common anions/inorganic elements
    inorganic_elements = {8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 
                         35: 'Br', 53: 'I', 7: 'N', 5: 'B', 14: 'Si'}
    
    for atom in atoms:
        atomic_num = atom.GetAtomicNum()
        formal_charge = atom.GetFormalCharge()
        
        # Check for metal ions
        if atomic_num in metals:
            has_metal = True
            metal_ions.add(metals[atomic_num])
            if formal_charge != 0:
                charged_atoms += 1
                
        # Count organic vs inorganic atoms
        elif atomic_num == 6:  # Carbon
            organic_atoms += 1
        elif atomic_num in inorganic_elements:
            inorganic_atoms += 1
            if formal_charge != 0:
                charged_atoms += 1
                
    # Check for water molecules
    water_pattern = Chem.MolFromSmiles("O")
    water_count = len(mol.GetSubstructMatches(water_pattern))
    
    # Rules for classification
    
    # Must have at least one metal
    if not has_metal:
        return False, "No metal atoms found"
    
    # Check for predominantly organic structure
    if organic_atoms > 8 and organic_atoms > inorganic_atoms:
        # Exception for simple organic salts (like magnesium stearate)
        if len(metal_ions) > 0 and charged_atoms >= 2:
            carboxylate = Chem.MolFromSmarts('C(=O)[O-]')
            if mol.HasSubstructMatch(carboxylate):
                return True, f"Metal organic salt containing {', '.join(metal_ions)}"
        return False, "Too many organic atoms for a mineral"
    
    # Should have some ionic character or coordination bonds
    if charged_atoms < 1 and water_count == 0:
        return False, "No ionic character or coordination bonds found"
    
    # Classify specific types of minerals
    if water_count > 0:
        return True, f"Hydrated mineral containing {', '.join(metal_ions)}"
    elif charged_atoms >= 2:
        return True, f"Ionic mineral containing {', '.join(metal_ions)}"
    else:
        return True, f"Mineral containing {', '.join(metal_ions)}"