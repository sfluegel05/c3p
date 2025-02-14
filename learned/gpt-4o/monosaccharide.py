"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones with three or more carbon
    atoms, potentially including cyclic forms like furanose or pyranose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has at least 3 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, "Too few carbon atoms for a monosaccharide"

    # Check for ring structures typical of monosaccharides (5 or 6 membered rings)
    rings = mol.GetRingInfo().AtomRings()
    has_valid_ring = any(len(ring) in {5, 6} for ring in rings)
    
    # Ensure multiple hydroxyl groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < min(c_count / 2, 3):  # Expect at least half carbons are hydroxylated
        return False, "Too few oxygen atoms for typical monosaccharide hydroxylation"

    # Look for potential carbonyl groups common in monosaccharide backbones
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)

    # Allow more complex structures, focusing on the anomalies that might arise in substituted monosaccharides
    # Examine potential functional groups typical in sugar derivatives only if bases are met
    if c_count > 12 or o_count > 8:
        # check for complex but relevant modifications like uronic acids, deoxy variants
        complex_patterns = [
            Chem.MolFromSmarts("[CX3;D2][OX2][CX3](=[OX1])"),
            Chem.MolFromSmarts("[CX4H2]"),  # Deoxy patterns 
        ]
        if not any(mol.HasSubstructMatch(pat) for pat in complex_patterns):
            return False, "Too complex for a simple monosaccharide without matching common derivations"

    # Confirm that the molecule is not overly complex with respect to known substituents
    molecular_formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    if "C" not in molecular_formula or "O" not in molecular_formula:
        return False, "Non-standard elemental composition for monosaccharides"

    # Classification reasoning
    if has_valid_ring or has_carbonyl:
        return True, "Structure has characteristic features of monosaccharides (polyhydroxy and cyclic or carbonyl)"
    
    return False, "Does not fit typical monosaccharide criteria"