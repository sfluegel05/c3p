"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: Xanthophylls (oxygenated carotenes), a subclass of carotenoids.
Heuristic: xanthophyll molecules are expected to have:
  - A long conjugated polyene chain (e.g. several adjacent conjugated double bonds),
  - A large hydrocarbon backbone (we require at least about 30 carbons),
  - At least one oxygenated functional group (a hydroxyl, carbonyl, or epoxide).
Note: This is a heuristic; false positives/negatives may still occur.
"""

from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids; they contain a long conjugated
    polyene system, a large carbon backbone (typically >30 carbons), and
    at least one oxygenated functional group (e.g. hydroxyl, carbonyl, or epoxide).

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

    # Count atoms to roughly assess size: carotenoids usually have a large carbon backbone (~C30+)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if carbon_count < 30:
        return False, "Too few carbons for a typical xanthophyll structure"
    if oxygen_count < 1:
        return False, "No oxygen atoms present, so not a xanthophyll"

    # Check for a long conjugated polyene chain.
    # First try a SMARTS pattern for at least 3 alternating double bonds.
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(polyene_pattern):
        # Fallback approach: count conjugated double bonds in the molecule.
        conjugated_dbl_bonds = 0
        for bond in mol.GetBonds():
            if (bond.GetBondTypeAsDouble() == 2.0) and bond.GetIsConjugated():
                conjugated_dbl_bonds += 1
        if conjugated_dbl_bonds < 5:
            return False, "No sufficiently long conjugated polyene chain detected"

    # Check for oxygenated functional groups.
    # Look for hydroxyl (-OH)
    hydroxyl = Chem.MolFromSmarts("[OX2H]")
    # Look for carbonyl (C=O)
    carbonyl = Chem.MolFromSmarts("[CX3]=O")
    # Look for epoxide: any three-membered ring with one oxygen and two carbons.
    epoxide = Chem.MolFromSmarts("[#6]-1-[#8]-[#6]-1")
    if not (mol.HasSubstructMatch(hydroxyl) or mol.HasSubstructMatch(carbonyl) or mol.HasSubstructMatch(epoxide)):
        return False, "No oxygenated functional groups detected"

    return True, "Molecule has a long conjugated polyene chain, large carbon backbone, and oxygen functionalities, consistent with xanthophyll structure."

# Example usage:
# test_smiles = "OC1CC(C(=C(C1)C)/C=C/C(/C)=C/C=C/C(/C)=C/C=C/C=C(/C=C/C=C(/C=C/C1=C(C)C[C@@H](O)CC1(C)C)"  # simplified example
# result, reason = is_xanthophyll(test_smiles)
# print(result, reason)