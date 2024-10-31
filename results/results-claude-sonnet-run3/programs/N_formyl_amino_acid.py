from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_formyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-formyl amino acid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an N-formyl amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for N-formyl amino acid pattern
    # [CH](=O)[NH][CH]C(=O)[OH]
    pattern = Chem.MolFromSmarts('[CH](=O)[NH][CH]C(=O)[OH]')
    
    if not mol.HasSubstructMatch(pattern):
        # Try alternate pattern for N-formyl amino acid
        # C(=O)[NH][CH]C(=O)[OH] 
        pattern2 = Chem.MolFromSmarts('C(=O)[NH][CH]C(=O)[OH]')
        if not mol.HasSubstructMatch(pattern2):
            return False, "No N-formyl amino acid pattern found"
    
    # Additional validation - check that formyl H is present
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            # Check if carbon is part of formyl group
            if atom.GetDegree() == 2:  # Formyl C should have 2 bonds
                neighbors = atom.GetNeighbors()
                if any(n.GetSymbol() == 'O' and n.GetDegree() == 1 for n in neighbors) and \
                   any(n.GetSymbol() == 'N' for n in neighbors):
                    return True, "Valid N-formyl amino acid structure found"
                    
    return False, "No valid formyl group found"
# Pr=1.0
# Recall=0.7142857142857143