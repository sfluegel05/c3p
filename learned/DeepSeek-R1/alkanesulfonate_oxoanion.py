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
    Must contain a sulfonate group (-SO3^-) attached to a carbon atom that is part of an alkane.

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

    # Define sulfonate oxoanion pattern: carbon (not in ring) connected to sulfur with three oxygens (two double bonds, one negative charge)
    # Pattern: [C;!R]-[S](=[O])(=[O])[O-]
    sulfonate_pattern = Chem.MolFromSmarts("[C;!R]-[S](=[O])(=[O])[O-]")
    
    # Check for sulfonate group presence
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "No sulfonate (-SO3^-) group attached to non-ring carbon found"

    # Verify sulfur atom has correct bonding and carbon is saturated
    sulfur_matches = mol.GetSubstructMatches(sulfonate_pattern)
    for match in sulfur_matches:
        carbon_idx = match[0]
        sulfur_idx = match[1]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        
        # Check carbon is saturated (all single bonds)
        if carbon_atom.GetTotalDegree() != carbon_atom.GetTotalNumHs() + sum(1 for bond in carbon_atom.GetBonds() if bond.GetBondType() == Chem.BondType.SINGLE):
            return False, "Sulfonate carbon is not fully saturated"
        
        # Check sulfur has exactly 4 bonds (two double to O, one single to O-, one single to C)
        if sulfur_atom.GetDegree() != 4:
            return False, "Sulfur in sulfonate group has incorrect bonding"

    return True, "Contains sulfonate group (-SO3^-) attached to alkane carbon"