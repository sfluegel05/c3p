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
    Determines if a molecule is a mucopolysaccharide (glycosaminoglycan).
    Criteria: Alternating uronic acid and glycosamine units in a polysaccharide structure,
    often sulfated. Uronic acid must be in a cyclic form (e.g., glucuronic acid).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect uronic acid (cyclic form with carboxylic acid)
    # Pattern for pyranose ring with COOH (like glucuronic acid)
    uronic_pattern = Chem.MolFromSmarts("[C;R][C](=O)[OH]")
    uronic_matches = mol.GetSubstructMatches(uronic_pattern)
    has_uronic = len(uronic_matches) > 0

    # Detect glycosamine (hexosamine in a ring with NH2/OH groups)
    # Pattern for amino sugar like glucosamine (ring with NH2 and adjacent O)
    glycosamine_pattern = Chem.MolFromSmarts("[N;H2,H1][C;R][C;R][O;R]")
    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)
    has_glycosamine = len(glycosamine_matches) > 0

    # Check sulfate esters (optional but common)
    sulfate_pattern = Chem.MolFromSmarts("[O]-S(=O)(=O)[O]")
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    # Glycosidic bonds (O connecting two cyclic carbons)
    glycosidic_pattern = Chem.MolFromSmarts("[C;R]@[O]@[C;R]")
    glycosidic_count = len(mol.GetSubstructMatches(glycosidic_pattern))

    # Check for polysaccharide characteristics (multiple rings)
    ring_info = mol.GetRingInfo()
    sugar_rings = sum(1 for ring in ring_info.AtomRings() if len(ring) in [5,6])  # pyranose/furanose
    reasons = []
    
    if not has_uronic:
        reasons.append("No uronic acid units detected")
    if not has_glycosamine:
        reasons.append("No glycosamine units detected")
    if glycosidic_count < 2:  # At least two glycosidic bonds for a chain
        reasons.append("Insufficient glycosidic bonds for polysaccharide")
    if sugar_rings < 2:
        reasons.append("Not enough sugar rings for polysaccharide")

    # Core criteria: presence of both units, sufficient glycosidic bonds, and rings
    if (has_uronic and has_glycosamine and glycosidic_count >=2 and sugar_rings >=2):
        reason = "Contains uronic acid and glycosamine units with glycosidic bonds in polysaccharide structure"
        if has_sulfate:
            reason += " and sulfate esters"
        return True, reason
    
    return False, "; ".join(reasons) if reasons else "Does not meet mucopolysaccharide criteria"