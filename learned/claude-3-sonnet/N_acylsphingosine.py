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

    # Core sphingosine pattern with more flexibility:
    # - CH2OH group
    # - NH-C(=O) amide connection
    # - Double bond in chain
    # - Two hydroxyls (one can be modified)
    sphingosine_core = Chem.MolFromSmarts("[CH2][OX2H]-[CH]([NX3H0,NX3H1]-[CX3]=[OX1])-[CH]([OX2H,OX2])-[CH]=[CH]")
    if not mol.HasSubstructMatch(sphingosine_core):
        return False, "Missing sphingosine core structure"

    # Verify amide group with connected carbon chain
    amide_pattern = Chem.MolFromSmarts("[NX3]([CX4])[CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No properly connected amide group found"

    # Count key elements
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    # Minimum requirements for elements
    if c_count < 16:  # Minimum for shortest ceramides
        return False, "Carbon count too low for N-acylsphingosine"
    if o_count < 3:  # 2 hydroxyls + 1 amide oxygen
        return False, "Insufficient oxygen count"
    if n_count < 1:  # Must have amide nitrogen
        return False, "No nitrogen found"
    
    # Verify presence of long carbon chains
    long_chain_pattern = Chem.MolFromSmarts("[CH2X4,CHX3][CH2X4,CHX3][CH2X4,CHX3][CH2X4,CHX3][CH2X4,CHX3][CH2X4,CHX3]")
    chain_matches = len(mol.GetSubstructMatches(long_chain_pattern))
    if chain_matches < 2:  # Need at least 2 long chains (sphingosine base and N-acyl)
        return False, "Missing required long carbon chains"

    # Check for disqualifying features
    
    # Exclude complex glycosylated derivatives
    glycosidic = Chem.MolFromSmarts("[OX2]([#6])[#6;R]")
    if mol.HasSubstructMatch(glycosidic):
        return False, "Contains glycosidic bonds"

    # Allow phosphate only in specific positions (for some derivatives)
    phosphate = Chem.MolFromSmarts("[P](=[O])([O,N])[O,N]")
    if mol.HasSubstructMatch(phosphate):
        # Only allow if connected to primary alcohol
        phosphate_ok = Chem.MolFromSmarts("[CH2][O][P](=[O])([O,N])[O,N]")
        if not mol.HasSubstructMatch(phosphate_ok):
            return False, "Contains incorrectly positioned phosphate group"

    # Molecular properties checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 1200:
        return False, "Molecular weight outside expected range"

    # Count rotatable bonds to ensure flexibility characteristic of ceramides
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for N-acylsphingosine structure"

    return True, "Contains characteristic N-acylsphingosine structure with proper backbone and N-acyl chain"