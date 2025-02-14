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

    # Pattern to identify an amino acid core: N-C-C(=O)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CH2,CH]C(=O)")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid functional group found"
    
    # Patterns to identify aromatic rings connected to the amino acid backbone or side chain
    aromatic_patterns = [
        Chem.MolFromSmarts("c1ccccc1"),  # Simple benzene ring
        Chem.MolFromSmarts("c1cc[nH]c1"),  # Pyrrole present in histidine
        Chem.MolFromSmarts("c1cc[n+](Cc2ccccc2)n1"),  # Imidazole ring
        Chem.MolFromSmarts("[cH]1[cH][cH][cH][cH]1"),  # Expanded generic aromatic
        Chem.MolFromSmarts("c1ccccc1C"),  # Aromatic chain branches any position
    ]
    
    for pattern in aromatic_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains amino acid functional group with an aromatic ring side chain"
    
    return False, "No aromatic ring attached as a side chain to an amino acid"