"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:37548 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    N-acylsphingosines are composed of sphingosine having a fatty acyl group 
    attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for key structural features:
    
    # 1. Amide group (N-C(=O)-)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"
    
    # 2. Two hydroxyl groups characteristic of sphingosine
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Missing required hydroxyl groups"
        
    # 3. Long carbon chain with double bond (sphingosine backbone)
    alkene_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX2]=[CX2]~[CX4,CX3]")
    if not mol.HasSubstructMatch(alkene_pattern):
        return False, "Missing characteristic sphingosine backbone with double bond"
    
    # 4. Check for primary alcohol (CH2-OH)
    primary_alcohol = Chem.MolFromSmarts("[CX4H2][OX2H]")
    if not mol.HasSubstructMatch(primary_alcohol):
        return False, "Missing primary alcohol group"

    # Count carbons and check molecular weight
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:  # Minimum for sphingosine (C18) plus small acyl group
        return False, "Carbon count too low for N-acylsphingosine"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for N-acylsphingosine"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Too few rotatable bonds for N-acylsphingosine structure"

    # Check for basic sphingosine structure with attached acyl group
    sphingosine_pattern = Chem.MolFromSmarts("[CX4H2][OX2H][CX4H]([NX3H])[CX4H]([OX2H])[CX3H]=[CX3H]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Missing characteristic sphingosine structure"

    return True, "Contains sphingosine backbone with N-acyl group and characteristic structural features"