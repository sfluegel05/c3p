"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenes, containing oxygen atoms in addition to the typical carotenoid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, not a xanthophyll"

    # Check for the presence of a long polyene chain (carotenoid backbone)
    # A simple way is to count the number of conjugated double bonds
    conjugated_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetIsConjugated())
    if conjugated_double_bonds < 8:  # Arbitrary threshold for a carotenoid-like structure
        return False, "Insufficient conjugated double bonds for a carotenoid backbone"

    # Check molecular weight (xanthophylls are typically >500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a xanthophyll"

    # Check for specific oxygen-containing functional groups (hydroxyl, epoxy, carbonyl)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    epoxy_pattern = Chem.MolFromSmarts("[OX2]C1CC1")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_epoxy = mol.HasSubstructMatch(epoxy_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)

    if not (has_hydroxyl or has_epoxy or has_carbonyl):
        return False, "No hydroxyl, epoxy, or carbonyl groups found"

    return True, "Contains a carotenoid backbone with oxygenated functional groups"