"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    
    A limonoid is a highly oxygenated triterpenoid with a prototypical structure either
    containing or derived from a precursor with a 4,4,8-trimethyl-17-furanylsteroid skeleton.

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
    furan_pattern = Chem.MolFromSmarts("c1cocc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring present"
    
    # Check for triterpenoid skeleton
    limonoid_core = Chem.MolFromSmarts("[C@]12[C@@](C)(CC[C@@]3([H])[C@@]4(C)C=CC(=O)OC(C)(C)[C@]4([H])C[C@@H]([C@@]3([H])C)C1=O)[C@@H](C[C@@]2([H])C)c1ccoc1"
    if not mol.HasSubstructMatch(limonoid_core):
        return False, "Missing limonoid core skeleton"
    
    # Check for high oxygenation
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygens < 5:
        return False, "Not highly oxygenated (fewer than 5 oxygens)"
    
    # Check for methyl groups at specific positions
    methyl_pattern = Chem.MolFromSmarts("[C@@]1(C)[C@@]2([H])[C@@]3([H])[C@]4([H])[C@@]5([H])[C@@]6([H])CC[C@@]7([H])[C@@]6(C)C=CC(=O)OC(C)(C)[C@]7([H])[C@@H]([C@@]5([H])[C@@]4(C)C(=O)[C@@]3([H])[C@@]2([H])C[C@@]1([H])C)C")
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Missing methyl groups at required positions"
    
    # Check for molecular weight and rotatable bonds (indicative of triterpenoid size)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for triterpenoid"
    
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if num_rotatable_bonds < 5:
        return False, "Too few rotatable bonds for triterpenoid"
    
    return True, "Meets the structural requirements of a limonoid"