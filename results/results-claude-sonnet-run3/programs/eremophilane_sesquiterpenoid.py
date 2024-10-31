from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_eremophilane_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is an eremophilane sesquiterpenoid.
    These are sesquiterpenoids (C15 compounds) with an eremophilane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an eremophilane sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check carbon count - should be 15 for sesquiterpenoids
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count < 15:
        return False, f"Not a sesquiterpenoid - has {carbon_count} carbons instead of at least 15"

    # The core eremophilane structure with flexible bond types
    # Matches decalin core with proper connectivity and gem-dimethyl group
    core_pattern = Chem.MolFromSmarts('[C,c]1[C,c][C,c][C,c]2[C,c][C,c][C,c][C,c][C,c]1[C,c]2(C)C')
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Missing characteristic eremophilane core structure"

    # Check for the presence of at least 2 methyl groups
    methyl_pattern = Chem.MolFromSmarts('[CH3]')
    if len(mol.GetSubstructMatches(methyl_pattern)) < 2:
        return False, "Missing characteristic methyl groups"

    # Check for ring system
    ring_info = mol.GetRingInfo()
    if not any(len(ring) == 6 for ring in ring_info.AtomRings()):
        return False, "Missing proper ring system"

    return True, "Eremophilane sesquiterpenoid with characteristic core structure and substitution pattern"
# Pr=None
# Recall=0.0