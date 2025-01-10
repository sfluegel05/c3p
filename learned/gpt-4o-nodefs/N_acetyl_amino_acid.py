"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group attached to the nitrogen of an amino acid.

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
    
    # Look for acetyl group ('CC(=O)N') attached directly to nitrogen
    acetyl_pattern = Chem.MolFromSmarts('CC(=O)N')
    acetyl_matches = mol.GetSubstructMatches(acetyl_pattern)
    if not acetyl_matches:
        return False, "No acetyl group attached to nitrogen found"
    
    # Look for carboxylic acid group (could be: 'C(=O)[O-]' or 'C(=O)O') as part of backbone
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)[O,#8]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    # Ensure the amino group (N) is part of the central amino acid chain, not elsewhere
    amino_acid_pattern = Chem.MolFromSmarts('[NX3][C;R0][C;!R0]')  # Free N attached to a non-ring C
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_matches:
        return False, "No amino acid backbone identified"
    
    # Check the attachment of acetyl group is at an amino acid nitrogen site
    for match in acetyl_matches:
        # Check if the nitrogen in the acetyl group is linked as part of the amino acid.
        acetyl_nitrogen_idx = [idx for idx in match if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7]
        for idx in acetyl_nitrogen_idx:
            if any(mol.GetAtomWithIdx(nei.GetIdx()).GetSmarts() == '[C;R0;$(C=[OX1,OX2])]C(=O)[O,#8]' for nei in mol.GetAtomWithIdx(idx).GetNeighbors()):
                return True, "Contains acetyl group attached to amino acid nitrogen with carboxylic acid group"

    return False, "Acetyl group not correctly bound to amino acid nitrogen"