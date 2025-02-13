"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    
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

    # Check for phosphocholine head group
    phosphocholine = Chem.MolFromSmarts("[O-]P(=O)(OCC[N+](C)(C)C)O")
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "Missing phosphocholine head group"

    # Check for sphingosine backbone with specific stereochemistry
    # [C@@H] and [C@H] ensure correct stereochemistry
    # \C=C\ ensures trans double bond
    sphingosine_pattern = Chem.MolFromSmarts("[CH2]CCCCCCCCCCCC/C=C/[C@@H](O)[C@H](COP)NC(=O)")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Missing characteristic sphingosine d18:1 backbone structure"

    # Count double bonds in molecule
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_count < 1:
        return False, "Missing required double bond"

    # Check for amide group
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatches(amide_pattern):
        return False, "Missing amide bond"

    # Check for hydroxyl group on sphingosine
    hydroxyl_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H](COP)")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing or incorrect hydroxyl group on sphingosine"

    # Verify N-acyl chain
    acyl_pattern = Chem.MolFromSmarts("NC(=O)CCCCC")  # At least 6 carbons
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Missing N-acyl chain"

    # Count key atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count != 2:  # One from amide, one from choline
        return False, "Must have exactly 2 nitrogen atoms"
    if p_count != 1:
        return False, "Must have exactly 1 phosphorus atom"
    if o_count < 6:  # Minimum 6 oxygens (phosphate, hydroxyl, amide, phosphocholine)
        return False, "Insufficient oxygen atoms"

    # Additional check for correct carbon chain length in sphingosine part
    sphingosine_chain = Chem.MolFromSmarts("CCCCCCCCCCCCC/C=C/[C@@H](O)[C@H]")
    if not mol.HasSubstructMatch(sphingosine_chain):
        return False, "Incorrect sphingosine chain length (must be 18 carbons)"

    return True, "Contains sphingosine d18:1 backbone with phosphocholine head group and N-acyl chain"