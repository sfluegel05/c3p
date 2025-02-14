"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid includes an amino acid structure with an aromatic ring attached as a side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern to identify an amino acid: N-C-C(=O)-O
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][C]([C])[CX3](=O)[OX2H]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid functional group found"
    
    # Aromatic side chain: aromatic group attached to the central or the side chain carbon
    aromatic_side_chain_patterns = [
        Chem.MolFromSmarts("[NX3][C]([C])C(=O)[O][C](a)"),  # attached directly to alpha carbon
        Chem.MolFromSmarts("[NX3][C]([C][CX4,CX3](a))C(=O)[O]")  # aromatic one step away
    ]
    
    for pattern in aromatic_side_chain_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains amino acid functional group with an aromatic ring side chain"
    
    return False, "No aromatic ring attached as a side chain to an amino acid"