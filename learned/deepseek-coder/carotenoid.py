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

    # Check for the presence of a long polyene chain (at least 7 conjugated double bonds)
    # More flexible pattern that allows for substitutions and variations
    polyene_pattern = Chem.MolFromSmarts("[C,c]=[C,c]~[C,c]=[C,c]~[C,c]=[C,c]~[C,c]=[C,c]~[C,c]=[C,c]~[C,c]=[C,c]")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No long polyene chain found (at least 7 conjugated double bonds required)"

    # Check molecular weight - carotenoids typically have a molecular weight around 400-700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 750:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for carotenoids (400-750 Da)"

    # Count carbons - carotenoids typically have around 30-50 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30 or c_count > 50:
        return False, f"Found {c_count} carbons, expected 30-50 for carotenoids"

    # Count double bonds - carotenoids typically have many double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bond_count < 7:
        return False, f"Found only {double_bond_count} double bonds, expected at least 7 for carotenoids"

    # Check for typical carotenoid end groups (optional)
    end_group_pattern = Chem.MolFromSmarts("[C,c]~[C,c]~[C,c]([C,c])~[C,c]~[C,c]")
    if not mol.HasSubstructMatch(end_group_pattern):
        return False, "No typical carotenoid end groups found"

    return True, "Contains a long polyene chain with conjugated double bonds and typical carotenoid characteristics"