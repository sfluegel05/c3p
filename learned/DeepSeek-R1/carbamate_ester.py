"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:23004 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is an ester of carbamic acid (structure: R-O-C(=O)-N<).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carbamate ester pattern: O connected to a carbon, then C(=O), then N with valence 3
    # The pattern [C]-O-C(=O)-[N;v3] ensures the O is part of an ester group
    carbamate_pattern = Chem.MolFromSmarts("[C][OX2]C(=O)[N;v3]")

    # Check for the presence of the carbamate group
    if mol.HasSubstructMatch(carbamate_pattern):
        # Additional check to exclude cases where N is part of a urea-like structure
        # Ensure the nitrogen is not connected to another carbonyl group
        for match in mol.GetSubstructMatches(carbamate_pattern):
            if len(match) < 4:
                continue  # Invalid match
            n_idx = match[3]
            n_atom = mol.GetAtomWithIdx(n_idx)
            # Check if nitrogen is connected to any carbonyl groups (other than the one in the carbamate)
            for neighbor in n_atom.GetNeighbors():
                if neighbor.GetIdx() == match[2]:  # Skip the carbonyl C in the carbamate
                    continue
                if neighbor.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in neighbor.GetBonds()):
                    # Check if the neighbor is a carbonyl carbon
                    for bond in neighbor.GetBonds():
                        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetBeginAtom().GetAtomicNum() == 8 or bond.GetEndAtom().GetAtomicNum() == 8:
                            return False, "Nitrogen connected to another carbonyl group (not a carbamate ester)"
        return True, "Contains carbamate ester group (O-C(=O)-N)"
    else:
        return False, "No carbamate ester group found"