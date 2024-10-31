from rdkit import Chem
from rdkit.Chem import AllChem

def is_16_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16-hydroxy steroid (steroid with a hydroxy group at position 16).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 16-hydroxy steroid, False otherwise 
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS patterns for steroid core with 16-OH
    # Pattern includes cyclopentanoperhydrophenanthrene core with OH at position 16
    # Multiple patterns to catch different variations and stereochemistry
    patterns = [
        # Basic steroid core with 16-OH (alpha)
        '[C]12[C][C][C]3[C]([C][C][C]4[C][C][C]([C@@H](O))[C][C]4[C]3)[C][C][C]2[C]1',
        # Basic steroid core with 16-OH (beta)
        '[C]12[C][C][C]3[C]([C][C][C]4[C][C][C]([C@H](O))[C][C]4[C]3)[C][C][C]2[C]1',
        # More general pattern for 16-OH steroid
        '[C]12[C][C][C]3[C]([C][C][C]4[C][C][C]([CH](O))[C][C]4[C]3)[C][C][C]2[C]1',
        # Pattern for prednisolone-like structures
        '[H][C]12[C]([C](O))[C]([C])(O)[C](=O)[C][C]1[C][C][C]1[C]2[C][C][C]2=CC(=O)[C]=[C][C]12',
        # Pattern for more complex steroids
        '[C]12[C][C][C]3[C]([C][C][C]4[C][C]([OH1])[C][C][C]4[C]3)[C][C][C]2[C]1'
    ]

    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            # Additional check for OH group
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 8:  # Oxygen
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 1:  # Hydrogen
                            return True, "16-hydroxy steroid identified"

    return False, "Does not contain steroid core with 16-hydroxy group"
# Pr=None
# Recall=0.0