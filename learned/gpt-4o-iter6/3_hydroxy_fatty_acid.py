"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxy functional group at the beta- or 3-position
    from the carboxylic acid and is characterized by a long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Locate the carboxylic acid group using a SMARTS pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    for match in carboxy_matches:
        # Match carbon atoms attached to the carboxylic carbon
        carboxylic_carbon = match[0]
        for atom in mol.GetAtomWithIdx(carboxylic_carbon).GetNeighbors():
            if atom.GetAtomicNum() == 6:  # Look for a carbon atom
                for next_atom in atom.GetNeighbors():
                    if next_atom.GetAtomicNum() == 6 and next_atom.GetIdx() != carboxylic_carbon:  # Second carbon
                        for final_atom in next_atom.GetNeighbors():
                            # Look for the hydroxy group attached to the third carbon
                            if (final_atom.GetAtomicNum() == 8  # Oxygen for hydroxyl group
                                and any(neigh.GetAtomicNum() == 1 for neigh in final_atom.GetNeighbors())):
                                return True, "Contains a hydroxy group at the 3-position"

    return False, "No hydroxy group at the 3-position found"