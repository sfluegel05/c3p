"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: polyunsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A PUFA must have:
    - A carboxylic acid group
    - Multiple carbon-carbon double bonds
    - A sufficiently long carbon chain

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PUFA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Count carbon-carbon double bonds
    # Exclude those in carboxylic group
    cc_double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(cc_double_bond_pattern))
    
    if double_bonds < 2:
        return False, f"Found only {double_bonds} C=C double bonds, need at least 2"

    # Count carbons (should be fatty acid-like)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:  # Most PUFAs have 18-22 carbons, but some can be shorter
        return False, f"Carbon chain too short ({c_count} carbons)"

    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 170:  # Approximate minimum weight for a PUFA
        return False, "Molecular weight too low"
    
    # Optional: Check for reasonable chain length using rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain appears too rigid for a PUFA"

    # Check for reasonable ratio of carbons to oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:  # Must have at least the carboxylic oxygens
        return False, "Too few oxygens"
    if o_count > c_count/2:  # Rough check to ensure it's primarily a hydrocarbon
        return False, "Too many oxygens relative to carbons"

    # If we get here, it's likely a PUFA
    return True, f"Contains carboxylic acid group and {double_bonds} C=C double bonds in a chain of {c_count} carbons"