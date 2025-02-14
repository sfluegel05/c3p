"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanethiol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thiol pattern: sulfur with one hydrogen attached to an sp3 carbon (alkyl carbon)
    thiol_pattern = Chem.MolFromSmarts("[#16H1][CX4;!$(C=*);!$(C~[#7,#8,#15,#16,#17,#35,#53])]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)

    if not thiol_matches:
        return False, "Does not contain an -SH group attached to an alkyl carbon"

    # Check each match to ensure the carbon is part of an alkyl group
    for match in thiol_matches:
        sulfur_idx, carbon_idx = match

        # Get the carbon atom attached to sulfur
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)

        # Check that carbon is only connected to carbons or hydrogens
        is_alkyl = True
        for neighbor in carbon_atom.GetNeighbors():
            atomic_num = neighbor.GetAtomicNum()
            if atomic_num not in [1, 6]:  # Hydrogen or Carbon
                is_alkyl = False
                break

        if not is_alkyl:
            continue  # Check next match
        else:
            return True, "Contains a sulfanyl group (-SH) attached to an alkyl group"

    return False, "Thiol group is not attached to an alkyl carbon"