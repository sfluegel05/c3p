"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: aldose sugars
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    Aldoses are polyhydroxy aldehydes or their cyclic hemiacetal forms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic size checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 250:  # Most simple aldoses are under 250 Da
        return False, "Molecule too large for simple aldose"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, "Too few carbons for an aldose"
    if c_count > 8:  # Most aldoses are C3-C7
        return False, "Too many carbons for simple aldose"
    if o_count < 3:
        return False, "Too few oxygens for an aldose"

    # Check for aromatic rings (aldoses shouldn't have these)
    if rdMolDescriptors.CalcNumAromaticRings(mol) > 0:
        return False, "Contains aromatic rings"

    # Pattern for aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")
    
    # More flexible hemiacetal pattern
    hemiacetal_pattern = Chem.MolFromSmarts("[O;R]-[CH1](-[O])-[C,O]")
    
    # Pattern for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    
    # Count hydroxyl groups
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, "Insufficient hydroxyl groups"

    # Check ring count
    rings = mol.GetRingInfo()
    if rings.NumRings() > 1:  # Simple aldoses have 0 or 1 ring
        return False, "Too many rings for simple aldose"
    
    # If cyclic, check ring size
    if rings.NumRings() == 1:
        ring_sizes = [len(r) for r in rings.AtomRings()]
        if not any(size in [5,6] for size in ring_sizes):
            return False, "Ring size not typical for aldose sugars"

    # Check for either aldehyde (open form) or hemiacetal (cyclic form)
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)
    
    if not (has_aldehyde or has_hemiacetal):
        return False, "No aldehyde or hemiacetal group found"

    # Pattern for carbon chain with hydroxyls
    polyhydroxy_chain = Chem.MolFromSmarts("[CH2,CH1,CH0]-[CH2,CH1,CH0]-[CH2,CH1,CH0]")
    if not mol.HasSubstructMatch(polyhydroxy_chain):
        return False, "No suitable carbon chain found"

    # More lenient ratio check for deoxy sugars
    ratio_oh_to_c = hydroxyl_matches / c_count
    if ratio_oh_to_c < 0.3:  # Lowered from 0.4 to catch deoxy sugars
        return False, "Too few hydroxyls relative to carbons"

    if has_aldehyde:
        return True, "Contains aldehyde group with multiple hydroxyls"
    else:
        return True, "Contains cyclic hemiacetal form with multiple hydroxyls"