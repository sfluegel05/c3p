"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: CHEBI:33877 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count nitrogen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # An amine must have at least 1 nitrogen atom
    if n_count == 0:
        return False, "No nitrogen atoms present"
    
    # Check if nitrogen is attached to hydrogen or carbon
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    is_amine = False
    for n_atom in nitrogen_atoms:
        neighbors = [mol.GetAtomWithIdx(x).GetAtomicNum() for x in n_atom.GetNeighbors()]
        if 1 in neighbors and 6 in neighbors:  # 1 = H, 6 = C
            is_amine = True
            break
    
    if is_amine:
        return True, "Contains nitrogen atom(s) attached to hydrogen and carbon"
    else:
        return False, "Nitrogen atom(s) not attached to both hydrogen and carbon"