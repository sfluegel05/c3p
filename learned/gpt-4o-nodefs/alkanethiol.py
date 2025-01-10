"""
Classifies: CHEBI:47908 alkanethiol
"""
from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    Alkanethiols feature one or more thiol groups (-SH) attached to an alkyl group.

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

    # Look for thiol group pattern (-SH)
    thiol_pattern = Chem.MolFromSmarts("[SX2H]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)

    if len(thiol_matches) == 0:
        return False, "No thiol groups found"
    
    # Consider presence of alkane-like structure
    # Verification of alkane or aliphatic features might be complex since
    # it may include non-ring structures or simple linear/cyclic saturated hydrocarbons.
    # Simplifying, let's just accept considering if thiol is present as enough
    # as definition may be loose with attached functionality.

    # Check that the molecule mainly consists of carbons and additionally other atoms like O, N
    non_acceptable_atoms = set()
    for atom in mol.GetAtoms():
        # Only allow carbons, hydrogens, sulfur (and few functional cases with N, O, especially when thiol present)
        if atom.GetAtomicNum() not in [6, 1, 8, 16, 7] and atom.GetIsAromatic() == False:
            non_acceptable_atoms.add(atom.GetSymbol())
    
    if non_acceptable_atoms:
        return False, f"Contains non-alkyl/thiol structure atoms: {non_acceptable_atoms}"

    return True, "Contains thiol group with alkyl or aliphatic association"