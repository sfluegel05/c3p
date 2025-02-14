"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a prenol phosphate resulting from the formal condensation
    of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

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

    # Check for phosphate or diphosphate group attached via oxygen
    phosphate_smarts = "[O]-P(=O)([O])[O]"
    diphosphate_smarts = "[O]-P(=O)([O])-O-P(=O)([O])[O]"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)

    if not phosphate_matches and not diphosphate_matches:
        return False, "No phosphate or diphosphate group found attached via oxygen"

    # Combine matches
    phosphate_atoms = set()
    for match in phosphate_matches + diphosphate_matches:
        phosphate_atoms.update(match)

    # Identify terminal phosphate group attached to terminal allylic carbon
    terminal_phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetIdx() in phosphate_atoms and atom.GetAtomicNum() == 8:
            # Check if this oxygen is connected to phosphate and carbon
            neighbors = atom.GetNeighbors()
            is_connected_to_phosphate = any(
                nbr.GetAtomicNum() == 15 and nbr.GetIdx() in phosphate_atoms for nbr in neighbors)
            is_connected_to_carbon = any(
                nbr.GetAtomicNum() == 6 for nbr in neighbors)
            if is_connected_to_phosphate and is_connected_to_carbon:
                carbon_atom = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 6][0]
                # Check if carbon is terminal (only connected to one heavy atom besides oxygen)
                carbon_neighbors = [nbr for nbr in carbon_atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != atom.GetIdx()]
                if len(carbon_neighbors) == 1:
                    # Check if carbon is allylic (adjacent to a double bond)
                    allylic = False
                    for nbr in carbon_neighbors:
                        bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
                        if bond.GetBondType() == rdchem.BondType.SINGLE and nbr.GetAtomicNum() == 6:
                            # Check if neighbor carbon has a double bond
                            for nnbr in nbr.GetNeighbors():
                                if nnbr.GetIdx() == carbon_atom.GetIdx():
                                    continue
                                bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nnbr.GetIdx())
                                if bond.GetBondType() == rdchem.BondType.DOUBLE and nnbr.GetAtomicNum() == 6:
                                    allylic = True
                                    break
                            if allylic:
                                break
                    if allylic:
                        terminal_phosphate_found = True
                        break
    if not terminal_phosphate_found:
        return False, "Phosphate group not attached to terminal allylic carbon"

    # Check for polyprenol backbone (multiple isoprene units)
    # Isoprene unit SMARTS allowing for cis/trans
    isoprene_smarts = "[CH2]=[CH]-[CH2]-[CH2]"
    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    # The molecule should have at least 2 isoprene units to be considered a polyprenol
    if len(isoprene_matches) < 2:
        return False, f"Insufficient number of isoprene units: found {len(isoprene_matches)} units"

    # Optionally, check the total number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 15:
        return False, f"Insufficient number of carbon atoms for a polyprenol backbone: found {num_carbons} carbons"

    return True, "Molecule is a polyprenol phosphate with a polyprenol backbone and terminal phosphate group attached to an allylic carbon"