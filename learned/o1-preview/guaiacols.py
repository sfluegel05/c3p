"""
Classifies: CHEBI:134251 guaiacols
"""
from rdkit import Chem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over atoms to find aromatic carbons with -OH attached
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
            # Check if this carbon has an -OH group attached
            has_OH = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                    # Found hydroxyl group
                    has_OH = True
                    break
            if has_OH:
                # Now check for adjacent aromatic carbons with methoxy group attached
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                        # Check if neighbor has a methoxy group attached
                        for nbr2 in neighbor.GetNeighbors():
                            if nbr2.GetAtomicNum() == 8 and nbr2.GetDegree() == 2:
                                # Oxygen with degree 2 may be methoxy
                                attached_methyl = False
                                for nbr3 in nbr2.GetNeighbors():
                                    if nbr3.GetAtomicNum() == 6 and nbr3.GetDegree() == 1 and nbr3.GetTotalNumHs() == 3:
                                        # Found methyl group attached to oxygen
                                        attached_methyl = True
                                        break
                                if attached_methyl:
                                    return True, "Contains phenol with methoxy group at ortho-position"
    return False, "Does not contain phenol with methoxy group at ortho-position"