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
    
    # 1. Amide group (N-C(=O)-) - more general pattern
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"
    
    # 2. Core sphingosine structure with flexibility for stereochemistry
    # [CH2OH]-[CH]-[NH]-[C]=O backbone with nearby OH and double bond
    sphingosine_core = Chem.MolFromSmarts("[CH2][OH].[CH]([NH])[CH]([OH])*.*=*")
    if not mol.HasSubstructMatch(sphingosine_core):
        return False, "Missing core sphingosine structure"

    # 3. Check for required functional groups
    
    # Primary alcohol (CH2-OH)
    primary_alcohol = Chem.MolFromSmarts("[CH2][OH]")
    if not mol.HasSubstructMatch(primary_alcohol):
        return False, "Missing primary alcohol group"
    
    # Secondary alcohol
    secondary_alcohol = Chem.MolFromSmarts("[CH]([#6])[OH]")
    if not mol.HasSubstructMatch(secondary_alcohol):
        return False, "Missing secondary alcohol group"
    
    # Double bond in chain
    alkene = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(alkene):
        return False, "Missing double bond"

    # Count carbons and check molecular weight
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:  # Minimum for basic N-acylsphingosine
        return False, "Carbon count too low for N-acylsphingosine"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for N-acylsphingosine"

    # Count nitrogens - should have exactly one (in the amide)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count != 1:
        return False, f"Should have exactly 1 nitrogen, found {n_count}"

    # Count oxygens - should have at least 3 (2 OH groups + 1 C=O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen count for N-acylsphingosine"

    # Verify long chain nature
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for N-acylsphingosine structure"

    return True, "Contains sphingosine backbone with N-acyl group and characteristic structural features"