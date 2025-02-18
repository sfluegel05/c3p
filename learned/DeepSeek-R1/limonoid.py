"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: CHEBI:55498 limonoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a prototypical structure containing a 4,4,8-trimethyl-17-furanylsteroid skeleton.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for furan-like fragment (5-membered ring with oxygen, allowing substitutions)
    furan_pattern = Chem.MolFromSmarts('[O;r5]1:[c,c,C;r5]:[c,c,C;r5]:[c,c,C;r5]:[c,c,C;r5]:1')
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan-like ring detected"
    
    # High oxygenation check (>=6 oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, f"Insufficient oxygen atoms ({o_count} < 6)"
    
    # Molecular weight typical of triterpenoids (modified limonoids often 450-1000 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 450:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"
    
    # Check for triterpenoid skeleton characteristics (carbon count 26-35 typical)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 26:
        return False, f"Carbon count too low ({c_count}) for triterpenoid"
    
    # Check for multiple oxygen-containing groups (esters, lactones, ketones, etc.)
    oxygen_groups = 0
    # Esters
    ester = Chem.MolFromSmarts('[OX2][CX3](=[OX1])')
    oxygen_groups += len(mol.GetSubstructMatches(ester))
    # Lactones (ester in ring)
    lactone = Chem.MolFromSmarts('[O;R][C;R](=O)')
    oxygen_groups += len(mol.GetSubstructMatches(lactone))
    # Ketones
    ketone = Chem.MolFromSmarts('[CX3](=O)[#6]')
    oxygen_groups += len(mol.GetSubstructMatches(ketone))
    # Ethers
    ether = Chem.MolFromSmarts('[OD2]([#6])[#6]')
    oxygen_groups += len(mol.GetSubstructMatches(ether))
    
    if oxygen_groups < 3:
        return False, f"Insufficient oxygen-containing groups ({oxygen_groups})"
    
    # Check ring system complexity (at least 4 rings)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, f"Not enough rings ({n_rings}) for complex triterpenoid"
    
    return True, "Contains furan ring, high oxygenation, and triterpenoid characteristics"