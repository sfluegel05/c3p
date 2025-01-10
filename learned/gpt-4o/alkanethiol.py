"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl or alkene group (terminus).

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
    
    # Broaden check to allow attachment to terminus carbon [sp3 include alkene endings]
    alkyl_attachment_pattern = Chem.MolFromSmarts("([#6,#7,#8][SX2H])")
    if not mol.HasSubstructMatch(alkyl_attachment_pattern):
        return False, "-SH group is not correctly attached to an alkyl or terminus alkene group"
    
    # Avoid complex peptides (like chains with many amino-acids)
    potential_protein_pattern = Chem.MolFromSmarts("C(=O)N")  # look for amide bond
    if mol.HasSubstructMatch(potential_protein_pattern):
        return False, "Structure too complex, likely proteinaceous"

    return True, "-SH group found and suitably attached"

# Example usage
smiles_examples = [
    "SCC(CCCC)CC", "CCS", "SC(CCCCC)C", "S\\C=C(\\CC)/C", "SCCCCCCCCCS", "SC(CCC)C"
]

for smiles in smiles_examples:
    result, reason = is_alkanethiol(smiles)
    print(f"SMILES: {smiles} -> {result}: {reason}")