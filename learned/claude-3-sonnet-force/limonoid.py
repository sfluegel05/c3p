"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: CHEBI:38708 limonoid 

A limonoid is any triterpenoid that is highly oxygenated and has a prototypical 
structure either containing or derived from a precursor with a 4,4,8-trimethyl-17-furanylsteroid skeleton.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.

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
    
    # Check for tetracyclic or rearranged tetracyclic skeleton
    tetracyclic_pattern = Chem.MolFromSmarts("[C&r5,r6,r7]1~[C&r5,r6,r7]2~[C&r5,r6,r7]3~[C&r5,r6,r7]4~[C&r5,r6,r7]1~[C&r5,r6,r7]2~[C&r5,r6,r7]3~[C&r5,r6,r7]4")
    if not mol.HasSubstructMatch(tetracyclic_pattern):
        return False, "No tetracyclic or rearranged tetracyclic skeleton found"
    
    # Check for presence of methyl groups
    methyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() == 3)
    if methyl_count < 3:
        return False, "Not enough methyl groups present"
    
    # Check for furan ring (or other heterocycles)
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    heterocycle_pattern = Chem.MolFromSmarts("[A&r5,r6]")
    if not mol.HasSubstructMatch(furan_pattern) and not mol.HasSubstructMatch(heterocycle_pattern):
        return False, "No furan ring or other heterocycle found"
    
    # Check oxygenation level
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_ratio = n_oxygens / n_carbons if n_carbons > 0 else 0
    if oxygen_ratio < 0.2:
        return False, "Not highly oxygenated enough for a limonoid"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 800:
        return False, "Molecular weight outside typical range for limonoids"
    
    return True, "Contains a tetracyclic or rearranged tetracyclic skeleton, methyl groups, a heterocycle, and high oxygenation level"