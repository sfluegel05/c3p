"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: CHEBI:17336 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Mucopolysaccharides are polysaccharides composed of alternating uronic acid and glycosamine units,
    often with sulfate ester groups.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule matches mucopolysaccharide characteristics
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for uronic acid components (carboxylic acid attached to a ring)
    uronic_acid_pattern = Chem.MolFromSmarts('[C;R](=O)[OX2H1]')  # Carboxylic acid on ring atom
    has_uronic = mol.HasSubstructMatch(uronic_acid_pattern)

    # Check for glycosamine components (amino group on a ring)
    glycosamine_pattern = Chem.MolFromSmarts('[N;R]')  # Any nitrogen in a ring
    has_glycosamine = mol.HasSubstructMatch(glycosamine_pattern)

    # Check for sulfate ester groups (O-SO3)
    sulfate_pattern = Chem.MolFromSmarts('[OX2]S(=O)(=O)[OX2]')
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    # Check for polysaccharide structure (multiple glycosidic bonds)
    glycosidic_pattern = Chem.MolFromSmarts('[r]O[r]')  # Oxygen connecting two rings
    glycosidic_bonds = len(mol.GetSubstructMatches(glycosidic_pattern))
    is_polymer = glycosidic_bonds >= 2  # At least two glycosidic linkages

    reasons = []
    if not has_uronic:
        reasons.append("No uronic acid units detected")
    if not has_glycosamine:
        reasons.append("No glycosamine units detected")
    if not is_polymer:
        reasons.append("Insufficient glycosidic bonds for polysaccharide structure")

    # Core requirements: alternating units implied by presence of both components in polymer
    if has_uronic and has_glycosamine and is_polymer:
        reason = "Contains alternating uronic acid and glycosamine units"
        if has_sulfate:
            reason += " with sulfate ester groups"
        return True, reason
    
    return False, "; ".join(reasons) if reasons else "Does not match mucopolysaccharide criteria"