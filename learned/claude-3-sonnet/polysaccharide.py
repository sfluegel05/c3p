"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: polysaccharide
A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
Should contain more than 10 monosaccharide residues.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple[bool, str]: (is_polysaccharide, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sugar ring pattern (6-membered ring with multiple OH groups)
    sugar_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1][CR1][OR1]1")
    sugar_rings = len(mol.GetSubstructMatches(sugar_pattern))
    
    if sugar_rings < 10:
        return False, f"Found only {sugar_rings} sugar rings, need at least 10 for polysaccharide"

    # Look for glycosidic linkages (C-O-C between rings)
    glycosidic_pattern = Chem.MolFromSmarts("[CR1]-[OR1]-[CR1]")
    glycosidic_bonds = len(mol.GetSubstructMatches(glycosidic_pattern))
    
    if glycosidic_bonds < 9:  # Need n-1 linkages for n rings
        return False, f"Found only {glycosidic_bonds} glycosidic linkages"

    # Count hydroxyl groups (characteristic of sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxyl_groups = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_groups < 20:  # Expect multiple OH groups per sugar
        return False, f"Too few hydroxyl groups ({hydroxyl_groups}) for polysaccharide"

    # Check molecular weight (should be large for polysaccharides)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 1500:  # Approximate minimum for 10+ sugar units
        return False, f"Molecular weight ({mol_weight:.1f}) too low for polysaccharide"

    # Count oxygen atoms (should be many in polysaccharides)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 30:  # Expect multiple O atoms per sugar unit
        return False, f"Too few oxygen atoms ({o_count}) for polysaccharide"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:  # Expect ~6 carbons per sugar unit
        return False, f"Too few carbon atoms ({c_count}) for polysaccharide"

    # Calculate ring fraction (should be high for polysaccharides)
    ring_fraction = len(mol.GetRingInfo().AtomRings()) / mol.GetNumAtoms()
    if ring_fraction < 0.2:  # At least 20% of atoms should be in rings
        return False, f"Ring fraction ({ring_fraction:.2f}) too low"

    return True, f"Contains {sugar_rings} sugar rings connected by {glycosidic_bonds} glycosidic linkages"