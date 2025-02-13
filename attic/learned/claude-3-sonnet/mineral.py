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
              46: 'Pd', 33: 'As', 3: 'Li', 4: 'Be', 37: 'Rb', 81: 'Tl'}
    
    # Common anions/inorganic elements
    inorganic_elements = {8: 'O', 9: 'F', 15: 'P', 16: 'S', 17: 'Cl', 
                         35: 'Br', 53: 'I', 7: 'N', 5: 'B', 14: 'Si'}
    
    # Track heteroatoms and their charges
    heteroatom_charges = []
    aromatic_atoms = 0
    
    for atom in atoms:
        atomic_num = atom.GetAtomicNum()
        formal_charge = atom.GetFormalCharge()
        
        if atom.GetIsAromatic():
            aromatic_atoms += 1
            
        # Check for metal ions
        if atomic_num in metals:
            has_metal = True
            metal_ions.add(metals[atomic_num])
            if formal_charge != 0:
                charged_atoms += 1
                heteroatom_charges.append(formal_charge)
                
        # Count organic vs inorganic atoms
        elif atomic_num == 6:  # Carbon
            if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) > 1:
                organic_atoms += 1
        elif atomic_num in inorganic_elements:
            inorganic_atoms += 1
            if formal_charge != 0:
                charged_atoms += 1
                heteroatom_charges.append(formal_charge)
                
    # Check for water molecules
    water_pattern = Chem.MolFromSmiles("O")
    water_count = len(mol.GetSubstructMatches(water_pattern))
    
    # Check for complex organic groups
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    
    # Exclusion rules
    if aromatic_atoms > 0:
        return False, "Contains aromatic rings - likely organic"
    
    if ring_count > 2:
        return False, "Too many rings for a mineral"
        
    if organic_atoms > 4 and not any(x in ["Mg", "Ca", "Ba", "Sr"] for x in metal_ions):
        return False, "Too many C-C bonds for a mineral"
        
    # Must have either a metal or metalloid (As, Sb)
    if not has_metal and not any(x in metal_ions for x in ['As', 'Sb']):
        return False, "No metal or metalloid atoms found"
    
    # Simple metal halides and oxides
    simple_inorganic_pattern = all(a.GetAtomicNum() in list(metals.keys()) + [9, 17, 35, 53, 8] for a in atoms)
    if simple_inorganic_pattern:
        return True, f"Simple inorganic mineral containing {', '.join(metal_ions)}"
    
    # Check charge balance for ionic compounds
    if charged_atoms > 0:
        if sum(heteroatom_charges) != 0 and water_count == 0:
            return False, "Unbalanced charges in structure"
    
    # Classification logic
    if water_count > 0:
        if organic_atoms > 2 and not any(x in ["Mg", "Ca", "Ba", "Sr"] for x in metal_ions):
            return False, "Hydrated organic compound"
        return True, f"Hydrated mineral containing {', '.join(metal_ions)}"
    
    elif charged_atoms >= 2:
        if organic_atoms > 2 and not any(x in ["Mg", "Ca", "Ba", "Sr"] for x in metal_ions):
            return False, "Organic salt"
        return True, f"Ionic mineral containing {', '.join(metal_ions)}"
    
    elif simple_inorganic_pattern:
        return True, f"Covalent mineral containing {', '.join(metal_ions)}"
        
    return False, "Structure not consistent with mineral characteristics"