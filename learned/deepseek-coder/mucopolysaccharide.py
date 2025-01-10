"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Mucopolysaccharides are polysaccharides composed of alternating units of uronic acids and glycosamines,
    often with sulfate ester groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More general uronic acid pattern (carboxyl group attached to a ring)
    uronic_acid_pattern = Chem.MolFromSmarts("[CX4;R][CX3](=[OX1])[OX2H0]")
    if not mol.HasSubstructMatch(uronic_acid_pattern):
        return False, "No uronic acid pattern found"

    # More general glycosamine pattern (amino group attached to a ring)
    glycosamine_pattern = Chem.MolFromSmarts("[CX4;R][NX3H0,H1,H2]")
    if not mol.HasSubstructMatch(glycosamine_pattern):
        return False, "No glycosamine pattern found"

    # Check for sulfate ester groups (common in mucopolysaccharides)
    sulfate_pattern = Chem.MolFromSmarts("[OX2][SX4](=[OX1])(=[OX1])[OX2]")
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    # Check molecular weight - mucopolysaccharides are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for mucopolysaccharide"

    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20:
        return False, "Too few carbons for mucopolysaccharide"
    if o_count < 10:
        return False, "Too few oxygens for mucopolysaccharide"
    if n_count < 1:
        return False, "Too few nitrogens for mucopolysaccharide"

    # Check for carbohydrate-like structure (multiple hydroxyl groups)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetTotalNumHs() > 0 and atom.GetAtomicNum() == 8)
    if hydroxyl_count < 3:
        return False, "Too few hydroxyl groups for carbohydrate structure"

    reason = "Contains uronic acids and glycosamines"
    if has_sulfate:
        reason += " with sulfate ester groups"
    return True, reason