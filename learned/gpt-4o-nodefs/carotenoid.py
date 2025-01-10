"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids typically have long conjugated polyene chains and/or cyclic end groups with
    functional groups such as hydroxyl, carbonyl, or epoxide groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Revised SMARTS patterns for carotenoid features
    long_conjugated_chain_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=[C]-[C]=[C](=[C]-[C]){3,}")  # more flexible on chain length
    cyclic_end_group_pattern = Chem.MolFromSmarts("[C]1=C(C(=C1)[C](=[C]-[C]){3,})")  # cyclic with extended conjugation

    # Check for a long conjugated chain
    has_conjugated_chain = mol.HasSubstructMatch(long_conjugated_chain_pattern)
    has_cyclic_end_group = mol.HasSubstructMatch(cyclic_end_group_pattern)

    # Check for functional groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    epoxide_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
    has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
    
    # Basic criteria: Must have a polyene chain and either functional groups or cyclic end groups
    if has_conjugated_chain and (has_cyclic_end_group or has_hydroxyl or has_carbonyl or has_epoxide):
        return True, "Contains long conjugated polyene chain and carotenoid-like features, including possible cyclic structures"

    # General framework of carotenoids: Atom count, bond count, molecular weight for complexity
    atom_count = mol.GetNumAtoms()
    bond_count = mol.GetNumBonds()
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    if atom_count > 40 and bond_count > 45 and mol_wt > 550:
        return True, "Large molecular framework and weight typical of carotenoids"

    return False, "Lacks typical carotenoid characteristics or size"

# Example usage (replace with actual example SMILES for testing)
example_smiles = "CC(\C=C\C=C(C)\C=C\C=C(C)\C=C\[C@@H](O)C(C)(C)C)=C/C=C/C=C(C)\C=C/C=C(C)/C=C/C1=C(C)C(O)[C@@H](O)CC1(C)C"
print(is_carotenoid(example_smiles))