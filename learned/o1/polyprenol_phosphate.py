"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a prenol phosphate resulting from the formal condensation
    of the terminal allylic hydroxy group of a polyprenol with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal phosphate group attached via oxygen
    phosphate_smarts = "[O]-P(=O)([O])[O]"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found attached via oxygen"

    # Check that the phosphate group is at the terminal position
    # Identify the oxygen atom connected to the phosphate group
    phosphate_atom_indices = set()
    for match in phosphate_matches:
        # The first atom in the SMARTS pattern is the oxygen connected to phosphate
        phosphate_oxygen_idx = match[0]
        phosphate_atom_indices.update(match)
        # Check if this oxygen is connected to only one non-hydrogen atom (terminal)
        oxygen_atom = mol.GetAtomWithIdx(phosphate_oxygen_idx)
        neighbors = [nbr for nbr in oxygen_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(neighbors) != 1:
            continue  # Not a terminal phosphate
        # Check if the oxygen is connected to a carbon atom
        carbon_atom = neighbors[0]
        if carbon_atom.GetAtomicNum() != 6:
            continue  # Not connected to carbon
        # Check if the carbon is allylic (adjacent to a double bond)
        is_allylic = False
        for nbr in carbon_atom.GetNeighbors():
            if nbr.GetIdx() == oxygen_atom.GetIdx():
                continue
            for bond in mol.GetBonds():
                if (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) in [(carbon_atom.GetIdx(), nbr.GetIdx()), (nbr.GetIdx(), carbon_atom.GetIdx())]:
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        is_allylic = True
                        break
            if is_allylic:
                break
        if not is_allylic:
            continue  # Phosphate not attached to allylic carbon
        break
    else:
        return False, "Phosphate group not attached to terminal allylic carbon"

    # Check for polyprenol backbone (repeating isoprene units)
    # Simplified approach: count the number of isoprene units by counting repeating patterns
    # Define an isoprene unit as: C=C-C-C
    isoprene_smarts = "C=C-C-C"
    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Insufficient number of isoprene units: found {len(isoprene_matches)} units"

    # Optionally, check the total number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 15:
        return False, f"Insufficient number of carbon atoms for a polyprenol backbone: found {num_carbons} carbons"

    return True, "Molecule is a polyprenol phosphate with a polyprenol backbone and terminal phosphate group attached to an allylic carbon"