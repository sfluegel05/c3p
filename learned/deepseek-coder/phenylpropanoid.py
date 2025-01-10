"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:26195 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is an organic aromatic compound with a structure based on a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one aromatic system (not just benzene)
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if not aromatic_atoms:
        return False, "No aromatic system found"

    # More flexible phenylpropane pattern matching
    # Allows for variations in the propyl chain and different aromatic systems
    phenylpropane_patterns = [
        # Standard patterns
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4]"),  # Standard phenylpropane
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]=[CX4]"),      # Double bond in chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4]"),       # Shorter chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][CX4]"),  # Longer chain
        
        # More flexible patterns
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][*]"),  # Substituted chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4]=[CX4][*]"),      # Substituted double bond
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][*]"),       # Substituted shorter chain
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][CX4][*]"),  # Substituted longer chain
        
        # Patterns with different aromatic systems
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][c]"),  # Connected to another aromatic
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][O]"),  # Connected to oxygen
        Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[CX4][CX4][CX4][N]"),  # Connected to nitrogen
    ]

    has_phenylpropane = any(mol.HasSubstructMatch(pattern) for pattern in phenylpropane_patterns)
    if not has_phenylpropane:
        return False, "No phenylpropane-like skeleton found"

    # Check molecular weight - phenylpropanoids typically >150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for phenylpropanoid"

    # If all checks pass, classify as phenylpropanoid
    return True, "Contains phenylpropane-like skeleton and meets molecular weight requirements"