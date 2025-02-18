"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:37550 sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin has a sphingoid base with an amide-linked fatty acid and
    a phosphocholine group attached to the terminal hydroxyl.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphocholine group (more flexible pattern)
    phosphocholine_pattern = Chem.MolFromSmarts("N(C)(C)CCOP(=O)(O)")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check for amide linkage to fatty acid (more flexible pattern)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide groups, need exactly 1"

    # Check for sphingoid base (more flexible pattern)
    sphingoid_pattern = Chem.MolFromSmarts("[C;H1,H2][C;H1,H2](O)[C;H1,H2][N;H1,H2]")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base pattern found"

    # Check molecular weight - sphingomyelins typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for sphingomyelin"

    # Count carbons, nitrogens and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 30:
        return False, "Too few carbons for sphingomyelin"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen (amide group)"
    if o_count < 5:
        return False, "Too few oxygens for sphingomyelin"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 15:
        return False, "Chains too short for sphingomyelin"

    return True, "Contains sphingoid base with amide-linked fatty acid and phosphocholine group"