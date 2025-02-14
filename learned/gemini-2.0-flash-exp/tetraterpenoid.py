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

    # Check for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O])([O])")
    if mol.HasSubstructMatch(phosphate_pattern):
       return False, "Contains phosphate group, likely not a tetraterpenoid"

    # Check for typical tetraterpenoid backbone patterns (more relaxed)
    #This backbone pattern tries to find linear chains with conjugated double bonds and methyl groups
    backbone_pattern = Chem.MolFromSmarts("C[C](C)=C[C](C)=C[C](C)=C[C](C)=C[C](C)=C[C](C)=C[C](C)=C[C](C)=C")
    if not mol.HasSubstructMatch(backbone_pattern):
       return False, "Does not contain typical tetraterpenoid backbone pattern"

    # Exclude molecules with glycoside groups
    glycoside_pattern = Chem.MolFromSmarts("OC[C]1[C](O)[C](O)[C](O)[C](O)[C](O)1")
    if mol.HasSubstructMatch(glycoside_pattern):
        return False, "Contains glycoside groups, likely not a tetraterpenoid"


    # Exclude quinones more stringently
    quinone_pattern1 = Chem.MolFromSmarts("C1(=O)C=CC=C(C=C1)=O") #conjugated 6-membered ring with two carbonyls
    quinone_pattern2 = Chem.MolFromSmarts("C1(=O)C(C)=CC=C(C=C1)=O") #conjugated 6-membered ring with two carbonyls and methyl
    if mol.HasSubstructMatch(quinone_pattern1) or mol.HasSubstructMatch(quinone_pattern2):
       return False, "Contains quinone ring, likely not a tetraterpenoid"


    # Exclude molecules with long aliphatic chains typical of fatty acids (C16+)
    long_chain_pattern1 = Chem.MolFromSmarts("CCCCCCCCCCCCCCCC") #simple C16 chain
    long_chain_pattern2 = Chem.MolFromSmarts("C(=O)OCCCCCCCCCCCCCCCC")#C16 with ester
    if mol.HasSubstructMatch(long_chain_pattern1) or mol.HasSubstructMatch(long_chain_pattern2):
        return False, "Contains long aliphatic chain, likely not a tetraterpenoid"

    # Exclude molecules with large ester groups typical of mycolic acids
    large_ester_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCCCCCCCCCCCCCCC")
    if mol.HasSubstructMatch(large_ester_pattern):
        return False, "Contains large ester group, likely not a tetraterpenoid"

    #Check for long chain ethers, these might be polyprenols or similar molecules
    ether_chain_pattern = Chem.MolFromSmarts("C-O-CCCCCCCCCCCCCCCC")
    if mol.HasSubstructMatch(ether_chain_pattern):
        return False, "Contains long ether chain, likely not a tetraterpenoid"


    # Count carbons (no limit)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:
        return False, f"Too few carbons to be a tetraterpenoid. {c_count}"

    return True, "Meets tetraterpenoid criteria: C40 isoprenoid structure, multiple double bonds."