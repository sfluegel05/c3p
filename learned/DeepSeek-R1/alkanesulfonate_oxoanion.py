"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
An alkanesulfonate in which the carbon at position 1 is attached to R (hydrogens, carbon chain, or other groups).
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    Must contain a sulfonate group (-SO3^-) attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sulfonate oxoanion pattern: carbon connected to sulfur with three oxygens (two double bonds, one negative charge)
    # Pattern: [C]-[S](=[O])(=[O])[O-]
    sulfonate_pattern = Chem.MolFromSmarts("[C]-[S](=[O])(=[O])[O-]")
    
    # Check for sulfonate group presence
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "No sulfonate (-SO3^-) group attached to carbon found"

    # Verify sulfur atom has correct bonding
    sulfur_matches = mol.GetSubstructMatches(sulfonate_pattern)
    for match in sulfur_matches:
        sulfur_idx = match[1]  # Sulfur is the second atom in the SMARTS pattern
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        
        # Check sulfur has exactly 4 bonds (single bond to C, two double to O, one single to O-)
        if sulfur_atom.GetDegree() != 4:
            return False, "Sulfur in sulfonate group has incorrect bonding"

    return True, "Contains sulfonate group (-SO3^-) attached to carbon atom"