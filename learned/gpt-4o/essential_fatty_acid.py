"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are a subset of polyunsaturated fatty acids required in the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a terminal carboxylic acid
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) < 1:
        return False, "No terminal carboxylic acid found"
    
    # Calculate carbon count to ensure it's a reasonable length for an essential fatty acid
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, f"Carbon chain is too short: only {carbon_count} carbon atoms"
    
    # Attempt to find multiple cis-configured double bonds for polyunsaturation
    cis_double_bond_pattern = Chem.MolFromSmarts("C/C=C\\C")
    double_bond_count = len(mol.GetSubstructMatches(cis_double_bond_pattern))
    
    # Relax the condition: A minimum of 3 cis double bonds might be too strict
    if double_bond_count < 2:
        return False, f"Insufficient cis double bonds: only {double_bond_count} found"
    
    # Exclude complex structures, esters, and functional groups not typical of free fatty acids
    # Prevent classification of conjugated or complex branched/cyclic derivatives
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Ester groups are not typical of free essential fatty acids"
    
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count > 0:
        return False, "Rings not typical of essential fatty acids"
    
    heteroatom_exceptions = any(atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1 and atom.GetAtomicNum() != 8 for atom in mol.GetAtoms())
    if heteroatom_exceptions:
        return False, "Presence of heteroatoms not common in essential fatty acids"
    
    return True, "Matches essential fatty acid pattern with polyunsaturation and a terminal carboxylic acid"