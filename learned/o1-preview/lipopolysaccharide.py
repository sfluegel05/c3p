"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: lipopolysaccharide
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide consists of saccharide units (including specific sugars like heptose and Kdo)
    attached to long-chain fatty acids via ester or amide linkages.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize molecule
    Chem.SanitizeMol(mol)

    # Detect sugar rings (rings with oxygen heteroatoms)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    sugar_rings = []
    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_size = len(ring)
        # Check if ring size is 5 to 7 (common for sugars)
        if ring_size >= 5 and ring_size <= 7:
            o_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            c_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
            # Sugars have one oxygen in the ring and the rest carbons
            if o_count == 1 and c_count == ring_size - 1:
                sugar_rings.append(ring)

    if len(sugar_rings) == 0:
        return False, "No sugar rings found"

    # Detect long aliphatic chains (length ≥ 12 carbons)
    # Looking for chains of at least 12 carbons
    fatty_acid_pattern = Chem.MolFromSmarts("[C;D2;H2][C;D2;H2][C;D2;H2][C;D2;H2][C;D2;H2][C;D2;H2][C;D2;H2][C;D2;H2][C;D2;H2][C;D2;H2][C;D2;H2][C;D2;H2]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) == 0:
        return False, "No long aliphatic chains found"

    # Detect ester or amide linkages
    ester_pattern = Chem.MolFromSmarts("[$([CX3](=O)[OX2H1]),$([CX3](=O)[OX1-])]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][CX4]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(ester_matches) == 0 and len(amide_matches) == 0:
        return False, "No ester or amide linkages found"

    # Check for molecular weight
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 1000:
        return False, "Molecular weight too low for lipopolysaccharide"

    return True, "Contains sugar rings linked to long-chain fatty acids via ester or amide bonds"

__metadata__ = {
    'chemical_class': {
        'name': 'lipopolysaccharide',
        'definition': 'Liposaccharide natural compounds consisting of a trisaccharide repeating unit (two heptose units and octulosonic acid) with oligosaccharide side chains and 3-hydroxytetradecanoic acid units (they are a major constituent of the cell walls of Gram-negative bacteria).'
    }
}