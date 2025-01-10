"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a prenol phosphate resulting from the formal condensation 
    of the terminal allylic hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.

    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for phosphate group connected via ester linkage
    phosphate_pattern = Chem.MolFromSmarts("O[P](=O)(O)O")
    if not phosphate_pattern:
        return False, "Invalid phosphate SMARTS pattern"

    # Search for phosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Define SMARTS pattern for isoprene unit: C=C-C-C
    isoprene_pattern = Chem.MolFromSmarts("C(=C)CC")
    if not isoprene_pattern:
        return False, "Invalid isoprene SMARTS pattern"

    # Find all isoprene units in the molecule
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    num_isoprene_units = len(isoprene_matches)

    if num_isoprene_units < 4:
        return False, f"Only {num_isoprene_units} isoprene units found; need at least 4 for polyprenol"

    # Check if phosphate group is connected to the polyprenol chain
    for phosphate_match in phosphate_matches:
        phosphate_o_atom_idx = phosphate_match[0]  # Oxygen atom connected to phosphate
        o_atom = mol.GetAtomWithIdx(phosphate_o_atom_idx)
        # Check neighbors of this oxygen atom
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                # Found carbon connected to the phosphate oxygen
                # Now traverse the carbon chain
                visited = set()
                to_visit = [neighbor.GetIdx()]
                chain_atoms = set()
                while to_visit:
                    current_idx = to_visit.pop()
                    if current_idx in visited:
                        continue
                    visited.add(current_idx)
                    current_atom = mol.GetAtomWithIdx(current_idx)
                    if current_atom.GetAtomicNum() != 6:
                        continue  # Only consider carbon atoms
                    chain_atoms.add(current_idx)
                    for nbr in current_atom.GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx not in visited and nbr.GetAtomicNum() == 6:
                            to_visit.append(nbr_idx)
                # Count the number of isoprene units within this chain
                isoprene_in_chain = 0
                for isoprene_match in isoprene_matches:
                    if all(idx in chain_atoms for idx in isoprene_match):
                        isoprene_in_chain += 1
                if isoprene_in_chain >= 4:
                    return True, f"Molecule is a polyprenol phosphate with {isoprene_in_chain} isoprene units"
    return False, "Phosphate group not connected to polyprenol chain with sufficient isoprene units"