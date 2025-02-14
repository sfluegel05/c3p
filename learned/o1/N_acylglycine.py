"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is an N-acyl-amino acid in which the amino acid specified is glycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define N-acylglycine pattern
    # The pattern represents R-C(=O)-N(H)-CH2-COOH
    n_acylglycine_pattern = Chem.MolFromSmarts('[CX3](=O)[NX3H][CH2][CX3](=O)[OX1H]')
    if not mol.HasSubstructMatch(n_acylglycine_pattern):
        return False, "No N-acylglycine substructure found"
    
    # Confirm that the backbone is glycine (no side chains on alpha carbon)
    alpha_carbon = Chem.MolFromSmarts('[NX3][CH2][CX3](=O)[OX1H]')
    matches = mol.GetSubstructMatches(alpha_carbon)
    if not matches:
        return False, "No glycine backbone found"
    for match in matches:
        alpha_c_idx = match[1]
        alpha_c_atom = mol.GetAtomWithIdx(alpha_c_idx)
        # Check if alpha carbon (CH2 group) is bonded to more than two hydrogens (shouldn't be)
        if len([nbr for nbr in alpha_c_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]) > 2:
            return False, "Alpha carbon has substituents, not glycine"
    
    return True, "Molecule is an N-acylglycine"