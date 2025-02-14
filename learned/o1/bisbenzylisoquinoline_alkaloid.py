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

    # Define isoquinoline ring SMARTS pattern
    isoquinoline_smarts = 'c1ccc2cccnc2c1'  # Isoquinoline ring
    isoquinoline_pattern = Chem.MolFromSmarts(isoquinoline_smarts)
    if isoquinoline_pattern is None:
        return False, "Invalid isoquinoline SMARTS pattern"

    # Find all isoquinoline rings in the molecule
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    if len(isoquinoline_matches) < 2:
        return False, f"Found {len(isoquinoline_matches)} isoquinoline units, need at least 2"

    # Get the atom indices for each isoquinoline unit
    isoquinoline_units = []
    for match in isoquinoline_matches:
        isoquinoline_units.append(set(match))

    # Check for bridges between isoquinoline units
    bridges_found = False
    for i in range(len(isoquinoline_units)):
        for j in range(i + 1, len(isoquinoline_units)):
            unit_i_atoms = isoquinoline_units[i]
            unit_j_atoms = isoquinoline_units[j]
            # Check for bonds connecting atoms from different units
            for bond in mol.GetBonds():
                atom1_idx = bond.GetBeginAtomIdx()
                atom2_idx = bond.GetEndAtomIdx()
                atom1_in_i = atom1_idx in unit_i_atoms
                atom2_in_j = atom2_idx in unit_j_atoms
                atom1_in_j = atom1_idx in unit_j_atoms
                atom2_in_i = atom2_idx in unit_i_atoms
                if (atom1_in_i and atom2_in_j) or (atom1_in_j and atom2_in_i):
                    # Atoms are from different isoquinoline units
                    atom1 = mol.GetAtomWithIdx(atom1_idx)
                    atom2 = mol.GetAtomWithIdx(atom2_idx)
                    # Ether bridge (oxygen connecting carbons from different units)
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and \
                       (atom1.GetAtomicNum() == 8 or atom2.GetAtomicNum() == 8):
                        bridges_found = True
                        return True, "Contains two isoquinoline units linked by an ether bridge"
                    # Carbon-carbon bridge
                    elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE and \
                         atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                        bridges_found = True
                        return True, "Contains two isoquinoline units linked by a carbon-carbon bridge"
                    # Methylenedioxy bridge (O-CH2-O bridge)
                    elif atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 8:
                        # Check if both oxygens are connected to the same carbon (methylenedioxy group)
                        carbons1 = [nbr for nbr in atom1.GetNeighbors() if nbr.GetAtomicNum() == 6]
                        carbons2 = [nbr for nbr in atom2.GetNeighbors() if nbr.GetAtomicNum() == 6]
                        common_carbons = set(carbons1).intersection(carbons2)
                        for carbon in common_carbons:
                            if carbon.GetDegree() == 3:
                                # Found a methylenedioxy bridge
                                bridges_found = True
                                return True, "Contains two isoquinoline units linked by a methylenedioxy bridge"
    if not bridges_found:
        return False, "No appropriate bridges found between isoquinoline units"
    else:
        return True, "Contains two isoquinoline units linked by bridges"