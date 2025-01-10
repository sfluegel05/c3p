"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the sulfanyl group (-SH)
    sulfanol_pattern = Chem.MolFromSmarts("[SX2H]") 
    if not mol.HasSubstructMatch(sulfanol_pattern):
        return False, "No -SH group found"
    
    # Ensure the -SH group is attached to an alkyl group (saturated carbon)
    alkyl_attachment_pattern = Chem.MolFromSmarts("[CX4][SX2H]") 
    if not mol.HasSubstructMatch(alkyl_attachment_pattern):
        return False, "-SH group is not attached to an alkyl group"

    return True, "-SH group found and attached to an alkyl group"

# Example usage
smiles_examples = [
    "SCC(CCCC)CC", "CCS", "SC(CCCCC)C", "S\\C=C(\\CC)/C", "SCCCCCCCCCS", "SC(CCC)C"
]

for smiles in smiles_examples:
    result, reason = is_alkanethiol(smiles)
    print(f"SMILES: {smiles} -> {result}: {reason}")