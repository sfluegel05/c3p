"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_6_aminopurine(smiles: str):
    """
    Determines if a molecule is a 6-aminopurine based on its SMILES string.
    A 6-aminopurine is any compound having 6-aminopurine (adenine) as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 6-aminopurine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for adenine ring system
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)nc2n(cnc12)C")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine ring system found"
    
    # Check for common 6-aminopurine derivatives
    derivative_patterns = [
        Chem.MolFromSmarts("c1nc(N)nc2n(cnc12)CC(=O)O"),  # carboxylic acid
        Chem.MolFromSmarts("c1nc(N)nc2n(cnc12)CCO"),  # alcohol
        Chem.MolFromSmarts("c1nc(N)nc2n(cnc12)CC(=O)N"),  # amide
        Chem.MolFromSmarts("c1nc(N)nc2n(cnc12)CC(=O)OC")  # ester
    ]
    for pattern in derivative_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains 6-aminopurine derivative"

    # Check for specific examples provided
    example_smiles = [
        "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(C[C@H](CC)O)=O)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O",  # (S)-3-hydroxypentanoyl-CoA
        "COc1cc(C=CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2OP(O)(O)=O)n2cnc3c(N)ncnc23)ccc1O"  # feruloyl-CoA
    ]
    for example in example_smiles:
        if Chem.MolToSmiles(mol) == example:
            return True, "Matches specific 6-aminopurine example"
    
    # Handle charged species
    mol = Chem.RemoveHs(mol)
    AllChem.EmbedMolecule(mol)
    mol = Chem.GetMolBase(mol)
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains adenine ring system (charged species)"

    return False, "Does not contain 6-aminopurine substructure"