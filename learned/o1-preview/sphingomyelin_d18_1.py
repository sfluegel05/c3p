"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is sphingomyelin d18:1 based on its SMILES string.
    A sphingomyelin d18:1 has a sphingosine backbone (d18:1) with an N-acyl fatty acid
    and a phosphocholine head group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Sphingosine backbone (d18:1 sphingoid base)
    sphingosine_smarts = "[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]=[#6]-[#6]-[#7]-[#6]-[#8]"
    sphingosine_pattern = Chem.MolFromSmarts(sphingosine_smarts)
    
    # N-acyl amide linkage
    amide_smarts = "[NX3][CX3](=O)[#6]"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)
    
    # Phosphocholine head group
    phosphocholine_smarts = "[O-]P(=O)(OCC[N+](C)(C)C)O"
    phosphocholine_pattern = Chem.MolFromSmarts(phosphocholine_smarts)
    
    # Check for sphingosine backbone
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone (d18:1) found"
    
    # Check for N-acyl amide linkage
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No N-acyl amide linkage found"
    
    # Check for phosphocholine head group
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine head group found"
    
    # Additional checks for accurate backbone length and double bond position
    # Count carbons in the sphingosine backbone
    sphingosine_matches = mol.GetSubstructMatch(sphingosine_pattern)
    sphingosine_atoms = [mol.GetAtomWithIdx(idx) for idx in sphingosine_matches]
    c_count = sum(1 for atom in sphingosine_atoms if atom.GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Sphingosine backbone has {c_count} carbons, expected 18"
    
    # Check for trans double bond at position 4
    double_bond_smarts = "[#6]-[#6]=[#6]-[#6]"
    double_bond_pattern = Chem.MolFromSmarts(double_bond_smarts)
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if not double_bond_matches:
        return False, "No trans double bond at position 4 in sphingosine backbone"
    
    return True, "Molecule is sphingomyelin d18:1 with correct sphingosine backbone, N-acyl amide linkage, and phosphocholine head group"