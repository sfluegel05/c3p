"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: CHEBI:17336 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide (glycosaminoglycan).
    Criteria: Polysaccharide with alternating uronic acid and glycosamine units,
    connected by glycosidic bonds, often sulfated.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved uronic acid detection: pyranose ring with COOH at C5
    uronic_pattern = Chem.MolFromSmarts(
        "[C;R6]1[C;R6][C;R6][C;R6][C;R6]([C;R6]1)(C(=O)O)"
    )
    uronic_count = len(mol.GetSubstructMatches(uronic_pattern))

    # Glycosamine: pyranose ring with NH2/N-acetyl at C2
    glycosamine_pattern = Chem.MolFromSmarts(
        "[C;R6]1[C;R6][N;H2,H1][C;R6][C;R6][C;R6]1"
    )
    glycosamine_count = len(mol.GetSubstructMatches(glycosamine_pattern))

    # Sulfate esters (optional but common)
    sulfate_pattern = Chem.MolFromSmarts("[O]S(=O)(=O)[O]")
    sulfate_count = len(mol.GetSubstructMatches(sulfate_pattern))

    # Glycosidic bonds (O connecting two anomeric carbons)
    glycosidic_pattern = Chem.MolFromSmarts("[C;R6]@[O]@[C;R6]")
    glycosidic_count = len(mol.GetSubstructMatches(glycosidic_pattern))

    # Check for polysaccharide characteristics
    ring_info = mol.GetRingInfo()
    sugar_rings = sum(1 for ring in ring_info.AtomRings() if len(ring) in [5,6])
    reasons = []

    # Molecular weight filter (polysaccharides are typically >1000 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 1000:
        reasons.append("Molecular weight too low for polysaccharide")

    if uronic_count < 1:
        reasons.append("No uronic acid units detected")
    if glycosamine_count < 1:
        reasons.append("No glycosamine units detected")
    if glycosidic_count < 3:  # Require longer chain
        reasons.append("Insufficient glycosidic bonds for polysaccharide")
    if sugar_rings < 4:  # At least 4 sugar units (alternating)
        reasons.append("Insufficient sugar rings for polysaccharide")

    # Check if both units are present and in sufficient quantity
    if uronic_count >=1 and glycosamine_count >=1 and glycosidic_count >=3 and sugar_rings >=4 and mol_wt >=1000:
        reason = "Alternating uronic acid and glycosamine units in polysaccharide structure"
        if sulfate_count > 0:
            reason += " with sulfate esters"
        return True, reason

    return False, "; ".join(reasons) if reasons else "Does not meet mucopolysaccharide criteria"