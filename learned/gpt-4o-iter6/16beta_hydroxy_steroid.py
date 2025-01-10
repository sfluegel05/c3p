"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
from rdkit import Chem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid has a hydroxyl group at position 16 of the steroid backbone
    with a beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for steroid backbone, flexible to slight structural deviations
    steroid_pattern = Chem.MolFromSmarts("[#6]1-[#6R2]-[#6R2]-[#6R2]-2-[#6R2]-[#6R2]-[#6R2]3-[#6R2]-[#6R2]-[#6R2]-[#6R2]-[#6R2]4-[#6R2]-[#6R2]-[#6R2]-[#6R2]-2-[#6R2]-1-3-4") 
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No recognizable steroid backbone found"

    # Simplified pattern that looks for the presence of an alcohol group on position 16, beta-configuration
    # Simplifying to check for a connected alcohol group (OH) to a carbon with at least one methyl group
    # Requires flexible positioning based on backbone assumptions
    position_16_beta_hydroxy = Chem.MolFromSmarts("[#6][C@H](O)[#6]([#6])C1CC(C)(C)C")
    if not mol.HasSubstructMatch(position_16_beta_hydroxy):
        return False, "No 16beta-hydroxy group found"
    
    return True, "Contains a steroid backbone with a 16beta-hydroxy group"

# Example test-case (a valid 16beta-hydroxy steroid example)
example_smiles = "C1=C2C(CC[C@]3([C@@]4(C[C@@H]([C@@H]([C@]4(CC[C@@]32[H])C)O)O)[H])[H])=CC(=C1)O"  # Assuming valid steroids as SMILES
result, reason = is_16beta_hydroxy_steroid(example_smiles)
print(f"Classification: {result}, Reason: {reason}")