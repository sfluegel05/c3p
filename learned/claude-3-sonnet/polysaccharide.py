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

    # Look for sugar ring patterns - allow for more variations
    sugar_pattern1 = Chem.MolFromSmarts("[C]1[C][C][C][C][O]1")  # Basic pyranose
    sugar_pattern2 = Chem.MolFromSmarts("[C]1[C][C][C][C]([O])[O]1")  # Alternative form
    
    sugar_rings = (len(mol.GetSubstructMatches(sugar_pattern1)) + 
                  len(mol.GetSubstructMatches(sugar_pattern2)))
    
    if sugar_rings < 10:
        return False, f"Found only {sugar_rings} sugar rings, need at least 10 for polysaccharide"

    # Look for glycosidic linkages (C-O-C between rings), allow for variations
    glycosidic_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
    glycosidic_bonds = len(mol.GetSubstructMatches(glycosidic_pattern))
    
    if glycosidic_bonds < 9:  # Need n-1 linkages for n rings
        return False, f"Found only {glycosidic_bonds} glycosidic linkages"

    # Count hydroxyl groups (characteristic of sugars)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")  # Allow for any OH
    hydroxyl_groups = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_groups < 15:  # Lowered threshold
        return False, f"Too few hydroxyl groups ({hydroxyl_groups}) for polysaccharide"

    # Check molecular weight (should be large for polysaccharides)
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 1000:  # Lowered threshold
        return False, f"Molecular weight ({mol_weight:.1f}) too low for polysaccharide"

    # Count oxygen atoms (should be many in polysaccharides)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 20:  # Lowered threshold
        return False, f"Too few oxygen atoms ({o_count}) for polysaccharide"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Lowered threshold
        return False, f"Too few carbon atoms ({c_count}) for polysaccharide"

    # Calculate ring fraction (should be present but can be lower)
    ring_fraction = len(mol.GetRingInfo().AtomRings()) / mol.GetNumAtoms()
    if ring_fraction < 0.05:  # Lowered threshold significantly
        return False, f"Ring fraction ({ring_fraction:.2f}) too low"

    # Additional check for carbohydrate nature - C:O ratio should be close to 1:1
    if not (0.8 <= c_count/o_count <= 1.5):
        return False, f"C:O ratio ({c_count/o_count:.2f}) not typical for polysaccharide"

    return True, f"Contains {sugar_rings} sugar rings connected by {glycosidic_bonds} glycosidic linkages"