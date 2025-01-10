"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    A carotenoid is a tetraterpenoid (C40) derived from psi,psi-carotene, with a long polyene chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a long polyene chain (at least 8 conjugated double bonds)
    # Use a more flexible pattern to account for substitutions and variations
    polyene_pattern = Chem.MolFromSmarts("[C,c]=[C,c]-[C,c]=[C,c]-[C,c]=[C,c]-[C,c]=[C,c]-[C,c]=[C,c]-[C,c]=[C,c]")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No long polyene chain found (at least 8 conjugated double bonds required)"

    # Check molecular weight - carotenoids typically have a molecular weight around 500-650 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450 or mol_wt > 700:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for carotenoids (450-700 Da)"

    # Count carbons - carotenoids typically have around 40 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30 or c_count > 50:
        return False, f"Found {c_count} carbons, expected around 40 for carotenoids"

    # Check for functional groups like hydroxyl, carbonyl, etc. (optional)
    functional_group_pattern = Chem.MolFromSmarts("[OH,OX1,C=O]")
    if not mol.HasSubstructMatch(functional_group_pattern):
        return False, "No functional groups found (some carotenoids have hydroxyl or carbonyl groups)"

    return True, "Contains a long polyene chain with conjugated double bonds, typical of carotenoids"