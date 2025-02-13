"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for furan ring
    furan_pattern = Chem.MolFromSmarts('c1ccoc1')
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Count carbons and oxygens - limonoids are highly oxygenated
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for a limonoid"
    if o_count < 4:
        return False, "Not enough oxygen atoms - limonoids are highly oxygenated"

    # Check molecular weight - should be substantial as these are complex molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, "Molecular weight too low for a limonoid"

    # Count rings - limonoids have multiple rings
    ring_info = mol.GetRingInfo()
    ring_count = ring_info.NumRings()
    if ring_count < 4:
        return False, "Too few rings for a limonoid structure"

    # Check for oxygen-containing functional groups
    ester_pattern = Chem.MolFromSmarts('[#6]-C(=O)-O-[#6]')
    ketone_pattern = Chem.MolFromSmarts('[#6]-C(=O)-[#6]')
    alcohol_pattern = Chem.MolFromSmarts('[#6]-[OH]')
    
    functional_groups = 0
    if mol.HasSubstructMatch(ester_pattern):
        functional_groups += len(mol.GetSubstructMatches(ester_pattern))
    if mol.HasSubstructMatch(ketone_pattern):
        functional_groups += len(mol.GetSubstructMatches(ketone_pattern))
    if mol.HasSubstructMatch(alcohol_pattern):
        functional_groups += len(mol.GetSubstructMatches(alcohol_pattern))
        
    if functional_groups < 2:
        return False, "Insufficient oxygen-containing functional groups"

    # Calculate number of sp3 carbons - limonoids have many
    sp3_c = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[CX4]')))
    if sp3_c < 10:
        return False, "Too few sp3 carbons for a limonoid skeleton"

    # Count number of methyl groups - limonoids typically have several
    methyl_pattern = Chem.MolFromSmarts('[CH3]')
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_count < 2:
        return False, "Too few methyl groups"

    return True, "Matches limonoid characteristics: contains furan ring, highly oxygenated, complex polycyclic structure with appropriate molecular weight and functional groups"