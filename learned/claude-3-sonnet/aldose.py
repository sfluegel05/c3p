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

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, "Too few carbons for an aldose"
    if o_count < 3:
        return False, "Too few oxygens for an aldose"

    # Pattern for aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CH1](=O)")
    
    # Pattern for hemiacetal in cyclic sugars (O-C-O where one O is in ring)
    hemiacetal_pattern = Chem.MolFromSmarts("[O;R]-[CH1]-[O]")
    
    # Pattern for hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    
    # Count hydroxyl groups
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_matches < 2:
        return False, "Insufficient hydroxyl groups"

    # Check for either aldehyde (open form) or hemiacetal (cyclic form)
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)
    
    if not (has_aldehyde or has_hemiacetal):
        return False, "No aldehyde or hemiacetal group found"

    # Pattern for carbon chain with hydroxyls
    polyhydroxy_chain = Chem.MolFromSmarts("[CH2,CH1,CH0]-[CH2,CH1,CH0]-[CH2,CH1,CH0]")
    if not mol.HasSubstructMatch(polyhydroxy_chain):
        return False, "No suitable carbon chain found"

    # Additional check for ring size if cyclic
    rings = mol.GetRingInfo()
    if rings.NumRings() > 0:
        ring_sizes = [len(r) for r in rings.AtomRings()]
        if not any(size in [5,6] for size in ring_sizes):
            return False, "Ring size not typical for aldose sugars"

    # Check that most carbons have hydroxyls (typical for sugars)
    ratio_oh_to_c = hydroxyl_matches / c_count
    if ratio_oh_to_c < 0.4:  # Allowing for some deoxy sugars
        return False, "Too few hydroxyls relative to carbons"

    if has_aldehyde:
        return True, "Contains aldehyde group with multiple hydroxyls"
    else:
        return True, "Contains cyclic hemiacetal form with multiple hydroxyls"