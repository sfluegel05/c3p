"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:36747 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sugar ring patterns
    sugar_ring_patterns = [Chem.MolFromSmarts(x) for x in ['OC1C(O)C(O)C(O)C1', 'OC1C(O)C(O)C(O)CC1']]
    sugar_rings = sum(mol.HasSubstructMatch(pat) for pat in sugar_ring_patterns)
    if sugar_rings < 4:
        return False, f"Found less than 4 sugar rings (found {sugar_rings})"
    
    # Look for glycosidic bonds (O-C-O)
    glycosidic_bond_pattern = Chem.MolFromSmarts('[OX2]C[OX2]')
    glycosidic_bonds = len(mol.GetSubstructMatches(glycosidic_bond_pattern))
    if glycosidic_bonds < 3:
        return False, f"Found less than 3 glycosidic bonds (found {glycosidic_bonds})"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 12 or c_count > 32:
        return False, "Carbon count outside expected range for tetrasaccharide"
    if o_count < 6 or o_count > 16:
        return False, "Oxygen count outside expected range for tetrasaccharide"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 1200:
        return False, "Molecular weight outside expected range for tetrasaccharide"
    
    return True, "Contains 4 sugar rings connected via glycosidic bonds"