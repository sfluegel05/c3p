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

    # Define benzylisoquinoline unit SMARTS pattern with variable substitutions
    benzylisoquinoline_smarts = '[#6]-c1ccc2ccnc(c2c1)-[*]'  # Simplified benzylisoquinoline pattern
    benzylisoquinoline_pattern = Chem.MolFromSmarts(benzylisoquinoline_smarts)
    if benzylisoquinoline_pattern is None:
        return False, "Invalid benzylisoquinoline SMARTS pattern"

    # Find all benzylisoquinoline units in the molecule
    matches = mol.GetSubstructMatches(benzylisoquinoline_pattern)
    if len(matches) < 2:
        return False, f"Found {len(matches)} benzylisoquinoline units, need at least 2"

    # Get the atom indices for each benzylisoquinoline unit
    benzylisoquinoline_units = []
    for match in matches:
        benzylisoquinoline_units.append(set(match))

    # Check for bridges between benzylisoquinoline units
    bridges_found = False
    for i in range(len(benzylisoquinoline_units)):
        for j in range(i + 1, len(benzylisoquinoline_units)):
            unit_i_atoms = benzylisoquinoline_units[i]
            unit_j_atoms = benzylisoquinoline_units[j]
            # Check for bonds connecting atoms from different units
            for bond in mol.GetBonds():
                atom1_idx = bond.GetBeginAtomIdx()
                atom2_idx = bond.GetEndAtomIdx()
                if (atom1_idx in unit_i_atoms and atom2_idx in unit_j_atoms) or \
                   (atom1_idx in unit_j_atoms and atom2_idx in unit_i_atoms):
                    # Check if the bond is an ether bridge or methylenedioxy group
                    atom1 = mol.GetAtomWithIdx(atom1_idx)
                    atom2 = mol.GetAtomWithIdx(atom2_idx)
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        if atom1.GetAtomicNum() == 8 or atom2.GetAtomicNum() == 8:
                            # Oxygen bridge (ether)
                            bridges_found = True
                            bridge_type = "ether bridge"
                            return True, f"Contains two benzylisoquinoline units linked by {bridge_type}"
                        elif atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                            # Carbon-carbon bridge
                            bridges_found = True
                            bridge_type = "carbon-carbon bridge"
                            return True, f"Contains two benzylisoquinoline units linked by {bridge_type}"
                    # Check for methylenedioxy group
                    elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE and \
                         (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or \
                         (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8):
                        # Possible methylenedioxy bridge
                        neighbor_atoms1 = [nb.GetAtomicNum() for nb in atom1.GetNeighbors()]
                        neighbor_atoms2 = [nb.GetAtomicNum() for nb in atom2.GetNeighbors()]
                        if 8 in neighbor_atoms1 and 8 in neighbor_atoms2:
                            bridges_found = True
                            bridge_type = "methylenedioxy bridge"
                            return True, f"Contains two benzylisoquinoline units linked by {bridge_type}"

    if not bridges_found:
        return False, "No appropriate bridges found between benzylisoquinoline units"
    else:
        return True, "Contains two benzylisoquinoline units linked by bridges"