"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is defined as a compound with repeated substructures (monomers), which can
    differ by chain length, composition, or branching.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is identified as a polymer-like structure, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for repetitive substructures
    # Example SMARTS for a simple aliphatic chain, common in small synthetic polymers
    repeat_unit_pattern = Chem.MolFromSmarts("*~*~*~*")
    repeat_unit_matches = mol.GetSubstructMatches(repeat_unit_pattern)
    
    if len(repeat_unit_matches) > 5:  # arbitrarily choosing an occurrence threshold
        return True, f"Contains repeating substructural patterns: found {len(repeat_unit_matches)}"

    # Analyze functional diversity (flexible for mixtures variation in functional groups)
    unique_atoms = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if len(unique_atoms) > 5:  # diversity in elemental composition
        return True, f"Contains diverse elemental composition: {len(unique_atoms)} different element types"

    # Check for large molecular size (common in polymers without known repeat units)
    if mol.GetNumAtoms() > 50:  # arbitrary size threshold
        return True, f"Large molecular size with {mol.GetNumAtoms()} atoms, possible polymer"

    return False, "No explicit polymer-like features detected"