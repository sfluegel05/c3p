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
    phosphate_smarts = "OP(O)(=O)[O-]"
    diphosphate_smarts = "OP(O)(=O)OP(O)(=O)[O-]"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)

    if not phosphate_matches and not diphosphate_matches:
        return False, "No phosphate or diphosphate group found"

    # Identify the phosphate group
    phosphate_atoms = set()
    for match in phosphate_matches + diphosphate_matches:
        phosphate_atoms.update(match)

    # Find terminal carbon connected to phosphate group
    terminal_phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Phosphorus atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen connected to phosphorus
                    for o_neighbor in neighbor.GetNeighbors():
                        if o_neighbor.GetAtomicNum() == 6 and o_neighbor.GetDegree() == 2:
                            # Potentially terminal carbon
                            terminal_carbon = o_neighbor
                            # Check if this carbon is at the end of a polyprenol chain
                            path = Chem.rdmolops.GetShortestPath(mol, terminal_carbon.GetIdx(), atom.GetIdx())
                            # Exclude short paths (should be longer for polyprenols)
                            if len(path) > 5:
                                terminal_phosphate_found = True
                                break
                    if terminal_phosphate_found:
                        break
            if terminal_phosphate_found:
                break
    if not terminal_phosphate_found:
        return False, "Phosphate group not attached to terminal carbon of polyprenol chain"

    # Check for polyprenol backbone (multiple isoprene units)
    # Isoprene units are 5-carbon units with alternating double bonds
    isoprene_unit = Chem.MolFromSmarts("C(=C)C-C=C")
    num_isoprene_units = len(mol.GetSubstructMatches(isoprene_unit))
    if num_isoprene_units < 2:
        return False, f"Insufficient number of isoprene units: found {num_isoprene_units} units"

    # Optionally, check the total number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 15:
        return False, f"Insufficient number of carbon atoms for a polyprenol backbone: found {num_carbons} carbons"

    return True, "Molecule is a polyprenol phosphate with a polyprenol backbone and terminal phosphate group attached to an allylic carbon"