"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: CHEBI:89833 sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    Sphingomyelin d18:1 has a sphingosine backbone (d18:1), a phosphocholine group, and an amide-linked fatty acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sphingosine backbone (d18:1) pattern
    # d18:1 sphingosine has 18 carbons with a double bond and a hydroxyl group
    # The pattern is more flexible to match different conformations
    sphingosine_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]=[CHX4][CH2X4][CHX4]([OHX2])[CH2X4][CHX4]([NHX3])[CX3](=[OX1])")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone (d18:1) found"

    # Check for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])([OX2]CC[N+](C)(C)C)")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Check for amide bond linking fatty acid chain
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Check for fatty acid chain (long carbon chain attached to amide)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 3:
        return False, f"Missing fatty acid chain, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - sphingomyelin d18:1 typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for sphingomyelin d18:1"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 30:
        return False, "Too few carbons for sphingomyelin d18:1"
    if o_count < 5:
        return False, "Too few oxygens for sphingomyelin d18:1"

    return True, "Contains sphingosine backbone (d18:1), phosphocholine group, and amide-linked fatty acid chain"