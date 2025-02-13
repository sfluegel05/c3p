"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_N_acylphytosphingosine, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide group (N-C(=O))
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Look for phytosphingosine core with more flexible pattern
    # Allow for variations in connectivity while maintaining core structure
    phyto_core = Chem.MolFromSmarts("[OX2H1,OX2H0][#6]-[#6]([NX3])-[#6]([OX2H1,OX2H0])-[#6]([OX2H1,OX2H0])-[#6]")
    if not mol.HasSubstructMatch(phyto_core):
        return False, "No phytosphingosine core structure found"

    # Count hydroxyl groups (both free and substituted)
    oh_pattern = Chem.MolFromSmarts("[OX2H1,OX2H0]-[#6]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 3:
        return False, f"Insufficient hydroxyl groups (found {oh_matches}, need ≥3)"

    # Check for fatty acyl chain - more flexible pattern
    # Look for carbon chain attached to amide
    fatty_chain = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No suitable fatty acyl chain found"

    # Basic size checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, f"Insufficient carbon atoms (found {c_count}, need ≥18)"

    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 300:
        return False, f"Molecular weight too low ({mol_weight:.1f} Da)"

    # Check for possible sugar modifications
    sugar_pattern = Chem.MolFromSmarts("[#6]1-[#8]-[#6]-[#6]-[#6]-[#6]1")
    has_sugar = mol.HasSubstructMatch(sugar_pattern)

    # Additional checks for chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain length too short for N-acylphytosphingosine"

    # Check nitrogen count
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen found"
    if n_count > 3 and not has_sugar:
        return False, f"Too many nitrogens for basic structure (found {n_count})"

    # Check for long chain characteristic of sphingolipids
    long_chain = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "Missing characteristic long carbon chain"

    return True, "Contains phytosphingosine backbone with N-acyl group and required hydroxyl groups"