from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_fluorinated_steroid(smiles: str):
    """
    Determines if a molecule is a fluorinated steroid.
    A fluorinated steroid is defined as a steroid substituted with one or more fluorine atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fluorinated steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Check for presence of fluorine atoms
    has_fluorine = False
    num_fluorines = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'F':
            has_fluorine = True
            num_fluorines += 1
            
    if not has_fluorine:
        return False, "No fluorine atoms present"

    # Get the basic steroid scaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    
    # Check for characteristic steroid ring system
    rings = scaffold.GetRingInfo()
    if rings.NumRings() < 4:
        return False, "Does not contain minimum 4 rings required for steroid core"
        
    # Check for characteristic 6-6-6-5 ring pattern of steroids
    ring_sizes = [len(ring) for ring in rings.AtomRings()]
    
    # Need at least one 5-membered ring and three 6-membered rings
    num_5_rings = sum(1 for size in ring_sizes if size == 5)
    num_6_rings = sum(1 for size in ring_sizes if size == 6)
    
    if num_5_rings < 1 or num_6_rings < 3:
        return False, "Does not have characteristic steroid 6-6-6-5 ring pattern"

    # If we get here, it's a fluorinated steroid
    return True, f"Fluorinated steroid with {num_fluorines} fluorine atom(s)"
# Pr=1.0
# Recall=1.0