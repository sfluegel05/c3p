from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrrolecarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a pyrrole carrying a single carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a pyrrolecarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for pyrrole substructure
    pyrrole_pattern = Chem.MolFromSmarts('[nH]1cccc1')
    if not mol.HasSubstructMatch(pyrrole_pattern):
        return False, "No pyrrole substructure found"

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    if len(carboxylic_matches) == 0:
        return False, "No carboxylic acid group found"
    elif len(carboxylic_matches) > 1:
        return False, "Multiple carboxylic acid groups found"

    # Check if carboxylic acid is attached to pyrrole ring
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    carboxylic_carbons = {match[0] for match in carboxylic_matches}
    
    for pyrrole_match in pyrrole_matches:
        pyrrole_atoms = set(pyrrole_match)
        # Check neighbors of pyrrole atoms
        for atom_idx in pyrrole_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            neighbor_ids = {n.GetIdx() for n in atom.GetNeighbors()}
            if neighbor_ids & carboxylic_carbons:
                position = None
                # Determine position of carboxylic acid on pyrrole
                pyrrole_n_idx = pyrrole_match[0]  # nitrogen is always first in the SMARTS match
                if atom_idx == pyrrole_match[1]:
                    position = 2
                elif atom_idx == pyrrole_match[2]:
                    position = 3
                elif atom_idx == pyrrole_match[3]:
                    position = 4
                elif atom_idx == pyrrole_match[4]:
                    position = 5
                return True, f"Pyrrole with carboxylic acid at position {position}"
    
    return False, "Carboxylic acid group not directly attached to pyrrole ring"
# Pr=1.0
# Recall=1.0