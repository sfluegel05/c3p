"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: CHEBI:55498 limonoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for furan ring (17-furanyl)
    furan = Chem.MolFromSmarts('o1cccc1')
    if not mol.HasSubstructMatch(furan):
        return False, "No furan ring detected"
    
    # Check oxygen count (highly oxygenated)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Only {o_count} oxygen atoms, expected >=4"
    
    # Check for ester groups (common in limonoids)
    ester = Chem.MolFromSmarts('[OX2][CX3](=[OX1])')
    ester_matches = mol.GetSubstructMatches(ester)
    if len(ester_matches) < 1:
        # Some might have lactones instead, check for lactone (ester in a ring)
        lactone = Chem.MolFromSmarts('[O;R][C;R](=O)')
        if not mol.HasSubstructMatch(lactone):
            return False, "No ester or lactone groups found"
    
    # Check molecular weight (triterpenoids are typically >400 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecule too light ({mol_wt:.1f} Da)"
    
    # Check for multiple methyl groups (4,4,8-trimethyl is part of the skeleton)
    methyl = Chem.MolFromSmarts('[CH3]')
    methyl_matches = len(mol.GetSubstructMatches(methyl))
    if methyl_matches < 3:
        return False, f"Only {methyl_matches} methyl groups, expected at least 3"
    
    # Basic triterpenoid check: ~30 carbons (may vary with substituents)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:
        return False, f"Only {c_count} carbons, unlikely triterpenoid"
    
    return True, "Contains furan ring, high oxygenation, and ester/lactone groups"