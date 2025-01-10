"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
"""

from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    An N-hydroxy-alpha-amino-acid is an amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all carboxyl groups (carbon with double bond to O and single bond to O-H)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # For each carboxyl group, check for alpha carbon and N-hydroxy group
    for match in carboxyl_matches:
        carboxyl_c_idx = match[0]  # Index of the carboxyl carbon
        carboxyl_c = mol.GetAtomWithIdx(carboxyl_c_idx)

        # Find alpha carbon (carbon connected to carboxyl carbon)
        alpha_c = None
        for neighbor in carboxyl_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                alpha_c = neighbor
                break
        if alpha_c is None:
            continue  # No alpha carbon found, skip to next carboxyl group

        # Check if alpha carbon is connected to nitrogen
        nitrogen_atom = None
        for neighbor in alpha_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:  # Nitrogen atom
                nitrogen_atom = neighbor
                break
        if nitrogen_atom is None:
            continue  # No nitrogen found, skip to next carboxyl group

        # Check if nitrogen has at least one hydroxy group attached
        has_hydroxy = False
        for nbr in nitrogen_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # Oxygen atom
                # Check if oxygen is connected to hydrogen (hydroxyl group)
                is_hydroxyl = False
                for o_nbr in nbr.GetNeighbors():
                    if o_nbr.GetAtomicNum() == 1:  # Hydrogen atom
                        is_hydroxyl = True
                        break
                if is_hydroxyl:
                    has_hydroxy = True
                    break
        if has_hydroxy:
            return True, "Contains N-hydroxy-alpha-amino-acid structure"

    return False, "Does not contain N-hydroxy-alpha-amino-acid structure"