"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Check molecular weight - sphingomyelins typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for sphingomyelin"

    # Check for phosphocholine head group
    phosphocholine = Chem.MolFromSmarts("[O-]P(=O)(OCC[N+](C)(C)C)O")
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "Missing phosphocholine head group"

    # Check for sphingosine backbone
    # Less strict pattern to match the core structure
    sphingosine_core = Chem.MolFromSmarts("[CH2]~[CH]~[CH](O)~[CH](CN)~[CH2]OP")
    if not mol.HasSubstructMatch(sphingosine_core):
        return False, "Missing characteristic sphingosine backbone structure"

    # Check for trans double bond in sphingosine
    trans_db = Chem.MolFromSmarts("C/C=C/C")
    if not mol.HasSubstructMatch(trans_db):
        return False, "Missing trans double bond in sphingosine"

    # Check for amide group
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if len(mol.GetSubstructMatches(amide_pattern)) != 1:
        return False, "Must have exactly one amide bond"

    # Check for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[CH](O)[CH](CN)")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing or incorrect hydroxyl group"

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

    # Count carbons to verify chain lengths
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:  # Minimum expected carbons for smallest sphingomyelin
        return False, "Insufficient carbon atoms for sphingomyelin"

    # Check for long carbon chains
    long_chain = Chem.MolFromSmarts("CCCCCCCC")  # At least 8 carbons in a row
    if len(mol.GetSubstructMatches(long_chain)) < 2:  # Need at least 2 long chains
        return False, "Missing required long carbon chains"

    return True, "Contains sphingosine d18:1 backbone with phosphocholine head group and N-acyl chain"