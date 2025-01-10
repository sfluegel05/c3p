"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI:21547 N-acetyl-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an N-acyl-amino acid that has acetyl as the acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core amino acid structure: NH2-CH-COOH
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4H][CX3](=[OX1])[OX2H,OX1H0-]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid structure found"

    # Look for the acetyl group (CH3-CO-) attached to the nitrogen
    acetyl_nitrogen_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])[NX3;H2,H1,H0][CX4H]")
    if not mol.HasSubstructMatch(acetyl_nitrogen_pattern):
        return False, "No acetyl group attached to nitrogen in amino acid structure"

    # Verify that the acetyl group is attached to the nitrogen of the amino acid
    # by checking if the nitrogen is part of both patterns
    nitrogen_matches = set()
    for match in mol.GetSubstructMatches(amino_acid_pattern):
        nitrogen_matches.add(match[0])  # Nitrogen is the first atom in the pattern
    
    acetyl_nitrogen_matches = set()
    for match in mol.GetSubstructMatches(acetyl_nitrogen_pattern):
        acetyl_nitrogen_matches.add(match[2])  # Nitrogen is the third atom in the pattern
    
    # Check if there's a nitrogen that's part of both patterns
    if not nitrogen_matches.intersection(acetyl_nitrogen_matches):
        return False, "Acetyl group not attached to the nitrogen of the amino acid"

    return True, "Contains an acetyl group attached to the nitrogen of an amino acid structure"