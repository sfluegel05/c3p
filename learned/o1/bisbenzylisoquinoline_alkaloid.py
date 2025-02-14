"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    A bisbenzylisoquinoline alkaloid consists of two benzylisoquinoline units linked by ether bridges,
    direct carbon-carbon bridges, or methylenedioxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define isoquinoline ring as SMARTS pattern (flexible matching)
    isoquinoline_smarts = 'c1ccc2ncccc2c1'  # Simplified isoquinoline pattern
    isoquinoline_pattern = Chem.MolFromSmarts(isoquinoline_smarts)

    if isoquinoline_pattern is None:
        return False, "Invalid isoquinoline SMARTS pattern"

    # Find all isoquinoline units in the molecule
    matches = mol.GetSubstructMatches(isoquinoline_pattern)
    if len(matches) < 2:
        return False, f"Found {len(matches)} isoquinoline units, need at least 2"

    # Get list of atom indices for each isoquinoline unit
    isoquinoline_units = []
    for match in matches:
        isoquinoline_units.append(set(match))

    # Check for bridges between isoquinoline units
    bridges_found = False
    for i in range(len(isoquinoline_units)):
        for j in range(i+1, len(isoquinoline_units)):
            unit_i_atoms = isoquinoline_units[i]
            unit_j_atoms = isoquinoline_units[j]
            # Find atoms that are connected between units
            for bond in mol.GetBonds():
                atom1_idx = bond.GetBeginAtomIdx()
                atom2_idx = bond.GetEndAtomIdx()
                if (atom1_idx in unit_i_atoms and atom2_idx in unit_j_atoms) or \
                   (atom1_idx in unit_j_atoms and atom2_idx in unit_i_atoms):
                    atom1 = mol.GetAtomWithIdx(atom1_idx)
                    atom2 = mol.GetAtomWithIdx(atom2_idx)
                    # Check if the bridge involves oxygen or carbon atoms
                    if atom1.GetAtomicNum() in [6, 8] and atom2.GetAtomicNum() in [6, 8]:
                        bridges_found = True
                        bridge_type = "oxygen" if 8 in [atom1.GetAtomicNum(), atom2.GetAtomicNum()] else "carbon"
                        return True, f"Contains two isoquinoline units linked by {bridge_type} bridge"
    if bridges_found:
        return True, "Contains two isoquinoline units linked by bridges"
    else:
        return False, "No bridges found between isoquinoline units"