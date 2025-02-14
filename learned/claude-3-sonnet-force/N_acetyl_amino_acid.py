"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI:38756 N-acetyl-amino acid
An N-acyl-amino acid that has acetyl as the acyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.

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
    
    # Look for acetyl group attached to nitrogen
    acetyl_pattern = Chem.MolFromSmarts("CC(=O)N")
    acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)
    if not acetyl_matches:
        return False, "No acetyl group attached to nitrogen"
    
    # Look for amino group
    amino_pattern = Chem.MolFromSmarts("[N;X3][C;X4]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if not amino_matches:
        return False, "No amino group found"
    
    # Check if the nitrogen attached to the acetyl group is also part of the amino group
    for acetyl_N in [match[1] for match in acetyl_matches]:
        for amino_N in [match[0] for match in amino_matches]:
            if acetyl_N == amino_N:
                break
        else:
            continue
        break
    else:
        return False, "Acetyl group not attached to amino nitrogen"
    
    # Look for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Check if the amino group is not part of a ring (to exclude proline)
    amino_atom = mol.GetAtomWithIdx(amino_N)
    if amino_atom.IsInRing():
        return False, "Amino group is part of a ring structure"
    
    return True, "Contains acetyl group attached to nitrogen of an amino acid"