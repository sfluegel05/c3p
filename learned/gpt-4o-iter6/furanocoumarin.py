"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a furan ring fused with a coumarin structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for furan-like pattern (5-membered oxygen-containing ring)
    furan_like_patterns = [
        Chem.MolFromSmarts("o1cccc1"),  # Typical furan ring
        Chem.MolFromSmarts("C1=COC=C1")  # Alternate depiction of furan
    ]
    furan_found = any(mol.HasSubstructMatch(f_pattern) for f_pattern in furan_like_patterns)
    if not furan_found:
        return False, "No furan-like ring found"
    
    # Look for coumarin-like structure pattern (various configurations)
    coumarin_like_patterns = [
        Chem.MolFromSmarts("O=c1ccc2c(c1)occ2"),  # Common coumarin pattern
        Chem.MolFromSmarts("O=c1cc2ccc(oc2c1)O")  # Alternate rings
    ]
    coumarin_found = any(mol.HasSubstructMatch(c_pattern) for c_pattern in coumarin_like_patterns)
    if not coumarin_found:
        return False, "No coumarin-like structure found"
    
    # Cross-check could involve examining the proximity connection of detected rings
    furan_atoms = [atom.GetIdx() for f_pattern in furan_like_patterns for atom in mol.GetSubstructMatches(f_pattern)]
    coumarin_atoms = [atom.GetIdx() for c_pattern in coumarin_like_patterns for atom in mol.GetSubstructMatches(c_pattern)]

    # Check if any atom from furan rings is neighboring any atom from coumarin
    for f_atom in furan_atoms:
        for c_atom in coumarin_atoms:
            if mol.GetBondBetweenAtoms(f_atom, c_atom):
                return True, "Furan ring and coumarin structure are connected, indicating fusion"

    return False, "Furan and coumarin structures not properly fused or no proper substructure found"