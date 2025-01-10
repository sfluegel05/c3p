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

    # Check for glycosidic bonds - exclude glycosylated ceramides
    glycosidic = Chem.MolFromSmarts("[OX2]([#6])[#6;R]")
    if mol.HasSubstructMatch(glycosidic):
        return False, "Contains glycosidic bonds - likely a glycosylated ceramide"

    # More specific sphingosine core pattern including the double bond position
    # [CH2OH]-[CH]-[NH]-C(=O) with OH and double bond in specific positions
    sphingosine_core = Chem.MolFromSmarts("[CH2][OH]-[CH]([NH][C]=[O])-[CH]([OH])-[CH]=[CH]")
    if not mol.HasSubstructMatch(sphingosine_core):
        return False, "Missing specific sphingosine core structure"

    # Check for amide group with proper context
    amide_pattern = Chem.MolFromSmarts("[NX3H]([CX4])[CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No properly connected amide group found"

    # Count carbons in longest chain from amide
    long_chain = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6][#6][#6][#6][#6][#6]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "N-acyl chain too short"

    # Basic element counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 16:
        return False, "Carbon count too low for N-acylsphingosine"
    if o_count < 3:
        return False, "Insufficient oxygen count"
    
    # Allow for some nitrogen-containing modifications if core structure is intact
    if n_count > 5:
        return False, "Too many nitrogens for N-acylsphingosine"

    # Check for phosphate groups that might indicate complex derivatives
    phosphate = Chem.MolFromSmarts("[P](=[O])([O,N])[O,N]")
    if mol.HasSubstructMatch(phosphate):
        return False, "Contains phosphate group - likely a complex derivative"

    # Verify aliphatic nature and chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for N-acylsphingosine structure"

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1000:
        return False, "Molecular weight outside expected range"

    return True, "Contains characteristic N-acylsphingosine structure with proper backbone and N-acyl chain"