"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are C40 isoprenoids, often with modifications like methyl loss/rearrangement.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for other elements
    other_elements = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 8])
    if other_elements > 3:
        return False, f"Too many other elements present, likely not tetraterpenoid: {other_elements} != (H, C, O) "

    # Check carbon count. Must be around 40.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 36 or c_count > 48:
       return False, f"Carbon count {c_count} is not within expected range (36-48)."
    
     # Check for typical tetraterpenoid backbone pattern
    backbone_pattern1 = Chem.MolFromSmarts("C=C-C=C-C=C-C=C-C=C-C=C-C=C-C=C-C=C") # Linear chain
    backbone_pattern2 = Chem.MolFromSmarts("C1CC(C)CC(C)C1") # Cyclic
    backbone_pattern3 = Chem.MolFromSmarts("C1=CC(C)=CC(C)=C1") #Cyclic with double bonds
    if not (mol.HasSubstructMatch(backbone_pattern1) or mol.HasSubstructMatch(backbone_pattern2) or mol.HasSubstructMatch(backbone_pattern3)):
        return False, "Does not contain typical tetraterpenoid backbone pattern"
    
    # Exclude quinones (check for conjugated 6-membered ring with two carbonyls)
    quinone_pattern = Chem.MolFromSmarts("C1(=O)C=CC(=O)C=C1")
    if mol.HasSubstructMatch(quinone_pattern):
       return False, "Contains quinone ring, likely not a tetraterpenoid"

    # Exclude molecules with long aliphatic chains typical of fatty acids (C16+)
    long_chain_pattern = Chem.MolFromSmarts("C[C](C)CCCCCCCCCCCCCCCC")
    if mol.HasSubstructMatch(long_chain_pattern):
        return False, "Contains long aliphatic chain, likely not a tetraterpenoid"

    # Exclude molecules with large ester groups typical of mycolic acids
    large_ester_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCCCCCCCCCCCCC")
    if mol.HasSubstructMatch(large_ester_pattern):
        return False, "Contains large ester group, likely not a tetraterpenoid"
    
    #Check for long chain ethers, these might be polyprenols or similar molecules
    ether_chain_pattern = Chem.MolFromSmarts("C-O-CCCCCCCCCCCCCCCC")
    if mol.HasSubstructMatch(ether_chain_pattern):
        return False, "Contains long ether chain, likely not a tetraterpenoid"
    
    
    # Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450:
        return False, f"Molecular weight too low: {mol_wt}. Requires at least 450"


    return True, "Meets tetraterpenoid criteria: C40 isoprenoid structure, multiple double bonds."