"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: CHEBI:38706 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclopenta[a]phenanthrene scaffold
    scaffold_pattern = Chem.MolFromSmarts("*1***2***3***4***5***6***7*1****6*****2*****3*****4*****5*****7*")
    if not mol.HasSubstructMatch(scaffold_pattern):
        return False, "No cyclopenta[a]phenanthrene scaffold found"

    # Check for methyl groups at C-10 and C-13
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    ring_info = mol.GetRingInfo()
    methyl_count = 0
    for match in methyl_matches:
        atom = mol.GetAtomWithIdx(match)
        if ring_info.IsBondInRingOfSize(atom.GetIdx(), 5):
            methyl_count += 1
    if methyl_count != 2:
        return False, f"Found {methyl_count} methyl groups in 5-membered rings, need exactly 2"

    # Check for optional alkyl group at C-17
    alkyl_pattern = Chem.MolFromSmarts("[C;D4]")
    alkyl_match = mol.GetSubstructMatches(alkyl_pattern)
    if len(alkyl_match) > 1:
        return False, "Found multiple alkyl groups, expected at most 1"

    # Check molecular weight (250-600 Da typical for steroids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:
        return False, "Molecular weight outside typical range for steroids"

    # Additional checks
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count < 3 or ring_count > 5:
        return False, "Number of rings outside typical range for steroids"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17 or c_count > 30:
        return False, "Number of carbon atoms outside typical range for steroids"

    return True, "Molecule contains cyclopenta[a]phenanthrene scaffold and meets other structural criteria for steroids"