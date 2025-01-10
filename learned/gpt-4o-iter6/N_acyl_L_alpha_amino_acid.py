"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    This class is defined by an L-alpha-amino acid backbone with an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # General pattern for L-alpha-amino acid backbone with stereochemistry specified
    l_alpha_amino_acid_pattern = Chem.MolFromSmarts('[C@H](N)([CX4])C(=O)O')  # Chiral carbon with amine and carboxyl connected
    if not mol.HasSubstructMatch(l_alpha_amino_acid_pattern):
        return (False, "No L-alpha-amino acid backbone found")

    # Pattern for nitrogen with acyl group
    n_acyl_pattern = Chem.MolFromSmarts('N([CX3](=O)[CX4,CX3])')  # N attached to carbonyl and possibly a chain
    n_acyl_matches = mol.GetSubstructMatches(n_acyl_pattern)

    if not n_acyl_matches:
        return (False, "No N-acyl group found attached to nitrogen")

    # Ensure N-acyl is on the alpha-amino nitrogen
    for match in n_acyl_matches:
        nitrogen_idx = match[0]
        atom = mol.GetAtomWithIdx(nitrogen_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.HasProp('_CIPCode'):
                if neighbor.GetProp('_CIPCode') in ['S', 'R']:  # chiral center
                    return (True, "Contains L-alpha-amino acid backbone with an N-acyl substituent correctly attached")

    return (False, "N-acyl group not correctly integrated into L-alpha-amino acid nitrogen")