"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: Xanthophylls (oxygenated carotenes), a subclass of carotenoids.
"""

from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids; they contain a long conjugated
    polyene system and also at least one oxygenated functional group
    (e.g. hydroxyl, carbonyl, or epoxide).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for the classification decision
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a heuristic SMARTS pattern for a conjugated polyene chain.
    # This pattern represents at least a chain with three alternating double bonds.
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No long conjugated polyene chain detected"

    # Define a SMARTS pattern for oxygenated functional groups:
    # This includes hydroxyl (-OH), carbonyl (C=O) and epoxide ([R3]O[R3]) groups.
    oxygen_pattern = Chem.MolFromSmarts("[$([OX2H]),$([CX3]=O),$([R3]O[R3])]")
    if not mol.HasSubstructMatch(oxygen_pattern):
        return False, "No oxygenated functional groups detected"

    # Optionally, one could add further checks such as the overall size of the conjugated system,
    # the location of the oxygen(s) relative to the polyene system, etc.
    return True, "Molecule has a conjugated polyene chain and oxygen functionalities, consistent with xanthophyll structure."

# Example usage (uncomment to test):
# test_smiles = "OC1CC(C(=C(C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C=C/C1=C(C)C[C@@H](O)CC1(C)C)"  # simplified
# result, reason = is_xanthophyll(test_smiles)
# print(result, reason)