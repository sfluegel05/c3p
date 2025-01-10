"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is a terpenoid derived from a monoterpene (C10 skeleton).
    The term includes compounds in which the C10 skeleton of the parent monoterpene
    has been rearranged or modified by the removal of one or more skeletal atoms
    (generally methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Monoterpenoids are derived from monoterpenes (C10 skeleton)
    # Allow for modifications, so require at least 7 carbons
    if c_count < 7:
        return False, f"Too few carbons ({c_count}) for a monoterpenoid"
    
    # Calculate the molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120 or mol_wt > 300:
        return False, f"Molecular weight ({mol_wt:.2f} Da) out of typical range for monoterpenoids (120-300 Da)"
    
    # Define SMARTS patterns for common monoterpene cores and rearranged forms
    monoterpene_patterns = [
        # Acyclic monoterpenes (e.g., myrcene, geraniol)
        "CC=C(C)CCC=C(C)C",          # General acyclic monoterpene pattern
        
        # Monocyclic monoterpenes (e.g., limonene, menthol)
        "C1CCC(CC1)C(C)C",           # Cyclohexane ring with isopropyl group
        "C1=CC(CC1)C(C)C",           # Cyclohexene ring with isopropyl group
        
        # Bicyclic monoterpenes (e.g., pinene, camphor)
        "C1([CH2])CC2CC1CC2C",       # Bicyclic ring system

        # Rearranged monoterpenes
        "C1CC2CCC1C2C",              # Bicyclo[2.2.1]heptane core
    ]
    
    # Check for presence of monoterpene cores
    core_found = False
    for pattern in monoterpene_patterns:
        smarts = Chem.MolFromSmarts(pattern)
        if smarts is None:
            continue
        if mol.HasSubstructMatch(smarts):
            core_found = True
            break
    if not core_found:
        return False, "No monoterpene core structures detected"

    # Monoterpenoids often contain oxygen-containing functional groups
    # Check for presence of hydroxyl, carbonyl, carboxyl, ester, ether groups
    functional_group_patterns = [
        "[OX2H]",                 # Hydroxyl group
        "[CX3]=[OX1]",            # Carbonyl group (aldehyde or ketone)
        "[CX3](=O)[OX2H1]",       # Carboxylic acid
        "[CX3](=O)[OX2][CX4]",    # Ester
        "[OX2][CX4]",             # Ether
        "[SX2H]",                 # Thiol group
        "[NX3][CX3]=[OX1]",       # Amide group
    ]
    fg_found = False
    for fg_smarts in functional_group_patterns:
        fg_pattern = Chem.MolFromSmarts(fg_smarts)
        if fg_pattern is None:
            continue
        if mol.HasSubstructMatch(fg_pattern):
            fg_found = True
            break
    if not fg_found:
        return False, "No typical monoterpenoid functional groups found"

    return True, "Molecule has characteristics of a monoterpenoid"