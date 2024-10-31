from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrahydroxyflavone(smiles: str):
    """
    Determines if a molecule is a tetrahydroxyflavone (any hydroxyflavone carrying four hydroxy substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydroxyflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic flavone core pattern (more specific pattern for flavone)
    flavone_pattern = Chem.MolFromSmarts('O=C1c2ccccc2OC(c2ccccc2)=C1')
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "No flavone core found"

    # Pattern for hydroxy groups
    oh_pattern = Chem.MolFromSmarts('[OH]')
    
    # Pattern for aromatic carbons with OH
    aromatic_oh_pattern = Chem.MolFromSmarts('c[OH]')
    
    # Get matches
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    aromatic_oh_matches = mol.GetSubstructMatches(aromatic_oh_pattern)
    
    # Get the atoms that are part of the flavone core
    core_match = mol.GetSubstructMatch(flavone_pattern)
    core_atoms = set(core_match)
    
    # Count hydroxy groups that are attached to the flavone system
    hydroxyl_count = 0
    for match in aromatic_oh_matches:
        carbon = match[0]  # The aromatic carbon attached to OH
        # Check if this carbon is part of or connected to the flavone core
        if carbon in core_atoms or any(neighbor.GetIdx() in core_atoms for neighbor in mol.GetAtomWithIdx(carbon).GetNeighbors()):
            hydroxyl_count += 1

    # Check the count
    if hydroxyl_count == 4:
        return True, "Found flavone with exactly 4 hydroxy groups"
    elif hydroxyl_count < 4:
        return False, f"Only {hydroxyl_count} hydroxy groups found, 4 required"
    else:
        return False, f"Found {hydroxyl_count} hydroxy groups (more than 4)"
# Pr=None
# Recall=0.0