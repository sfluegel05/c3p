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
    
    # Check for alkyl group attachment (any carbon attached to -SH specifically)
    # Alkanethiols can have unsaturations or other substituents in the molecule
    attachment_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2][SX2H]")  # Carbon followed by sulfur
    if not mol.HasSubstructMatch(attachment_pattern):
        return False, "-SH group is not suitably attached to a carbon group"
    
    return True, "-SH group found and suitably attached to an alkyl group"

# Example usage
smiles_examples = [
    "SCC(CCCC)CC", "CCS", "SC(CCCCC)C", "S\\C=C(\\CC)/C", "SCCCCCCCCCS", "SC(CCC)C"
]

for smiles in smiles_examples:
    result, reason = is_alkanethiol(smiles)
    print(f"SMILES: {smiles} -> {result}: {reason}")