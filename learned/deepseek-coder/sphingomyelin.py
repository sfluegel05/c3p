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

    # More specific phosphocholine pattern that includes the ester linkage
    phosphocholine_pattern = Chem.MolFromSmarts("[N+](C)(C)CCOP(=O)([O-])O")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check for amide linkage to fatty acid
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide groups, need exactly 1"

    # More specific sphingoid base pattern that includes the characteristic structure
    sphingoid_pattern = Chem.MolFromSmarts("[C;H1,H2][C;H1,H2](O)[C;H1,H2][N;H1,H2]")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base pattern found"

    # Check molecular weight - sphingomyelins typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for sphingomyelin"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 25:
        return False, "Too few carbons for sphingomyelin"
    if o_count < 5:
        return False, "Too few oxygens for sphingomyelin"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short for sphingomyelin"

    # Additional check for the characteristic sphingomyelin structure
    # Should have at least one long carbon chain (>= 12 carbons)
    long_chain_pattern = Chem.MolFromSmarts("[C;H2][C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]~[C;H2]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain found"

    return True, "Contains sphingoid base with amide-linked fatty acid and phosphocholine group"