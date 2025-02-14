"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group on the carbon atom adjacent to (alpha to) the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Patterns for alpha-amino acid structure
    patterns = [
        Chem.MolFromSmarts("[NX3H2,NX3H1]-[CX4H1]-[CX3](=O)-[OX1-,OX2H1]"),    # Simple linear alpha-amino acid
        Chem.MolFromSmarts("[NX3H2,NX3H1]-[C@H]-[CX3](=O)-[OX1-,OX2H1]"),     # Chiral alpha carbon
        Chem.MolFromSmarts("[NX3H2,NX3H1]-[C@@H]-[CX3](=O)-[OX1-,OX2H1]"),    # Reverse chirality alpha carbon
        Chem.MolFromSmarts("[NX3H2,NX3H1]-*-[C](=O)-[OX1-,OX2H1]"),           # Generic pattern accommodating variations
        Chem.MolFromSmarts("[NX3,NX4]-[CX4]-[CX3]=O")                        # Rarer formal representations
    ]
    
    # Check the molecule against all patterns
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains alpha-amino acid structure"

    return False, "No alpha-amino acid pattern found"