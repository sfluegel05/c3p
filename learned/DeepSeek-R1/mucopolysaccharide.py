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
    Mucopolysaccharides are polysaccharides with alternating uronic acid and glycosamine units,
    often sulfated.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Uronic acid: carboxylic acid group (COOH)
    uronic_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    has_uronic = mol.HasSubstructMatch(uronic_pattern)

    # Glycosamine: amino group attached to carbon adjacent to oxygen (sugar-like)
    glycosamine_pattern = Chem.MolFromSmarts('[N;H2,H1][C][O]')
    has_glycosamine = mol.HasSubstructMatch(glycosamine_pattern)

    # Sulfate ester (O-SO3)
    sulfate_pattern = Chem.MolFromSmarts('[OX2]S(=O)(=O)[OX2]')
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    # Glycosidic bonds (O connecting two carbons with multiple oxygen attachments)
    glycosidic_pattern = Chem.MolFromSmarts('[C;!$(C=O)]O[C;!$(C=O)]')
    glycosidic_count = len(mol.GetSubstructMatches(glycosidic_pattern))

    reasons = []
    if not has_uronic:
        reasons.append("No uronic acid units detected")
    if not has_glycosamine:
        reasons.append("No glycosamine units detected")
    if glycosidic_count < 1:
        reasons.append("Insufficient glycosidic bonds")

    # Core requirements: both components and at least one glycosidic bond
    if has_uronic and has_glycosamine and glycosidic_count >= 1:
        reason = "Contains uronic acid and glycosamine units with glycosidic bonds"
        if has_sulfate:
            reason += " and sulfate esters"
        return True, reason
    
    return False, "; ".join(reasons) if reasons else "Does not meet mucopolysaccharide criteria"