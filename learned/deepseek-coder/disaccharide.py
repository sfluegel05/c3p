"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:36218 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound in which two monosaccharides are joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least two sugar-like rings
    # A sugar ring typically has 5 or 6 atoms, with multiple hydroxyl groups
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "No rings found, need at least one for a disaccharide"

    # Check for glycosidic bond (oxygen linking two sugar-like rings)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[O;X2][C;X4][C;X4][O;X2]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond found between rings"

    # Count hydroxyl groups (OH) to ensure it's a sugar
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if hydroxyl_count < 4:
        return False, "Too few hydroxyl groups for a disaccharide"

    # Check molecular weight - disaccharides typically have a molecular weight between 200 and 600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside typical range for disaccharides"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8 or c_count > 24:
        return False, "Carbon count outside typical range for disaccharides"
    if o_count < 6 or o_count > 14:
        return False, "Oxygen count outside typical range for disaccharides"

    return True, "Contains two monosaccharide units joined by a glycosidic bond"