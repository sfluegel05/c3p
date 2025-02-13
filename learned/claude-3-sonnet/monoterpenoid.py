"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:36925 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    Monoterpenoids are terpenoids derived from a monoterpene, including compounds where
    the C10 skeleton of the parent monoterpene has been rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check atom count range (8-25 atoms)
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 8 or n_atoms > 25:
        return False, f"Atom count ({n_atoms}) outside typical range for monoterpenoids"
    
    # Look for common monoterpenoid substructures and functional groups
    substructure_patterns = (
        Chem.MolFromSmarts("[C&r6,C&r5]"),  # Rearranged/modified monoterpene skeleton
        Chem.MolFromSmarts("[OX1]"),  # Alcohol
        Chem.MolFromSmarts("[OX2]"),  # Ether
        Chem.MolFromSmarts("[C&r6,C&r5][OX2]"),  # Epoxide
        Chem.MolFromSmarts("[CX3](=[OX1])"),  # Ketone
        Chem.MolFromSmarts("[C&r6,C&r5]=[C&r6,C&r5]"),  # Exocyclic double bond
        Chem.MolFromSmarts("[C&r6,C&r5]1[C&r6,C&r5][C&r6,C&r5][C&r6,C&r5]1"),  # Cyclobutane ring
        Chem.MolFromSmarts("[C&r6,C&r5]1[C&r6,C&r5][C&r6,C&r5][C&r6,C&r5][C&r6,C&r5]1"),  # Cyclopentane ring
        Chem.MolFromSmarts("[C&r6,C&r5]1[C&r6,C&r5][C&r6,C&r5][C&r6,C&r5][C&r6,C&r5][C&r6,C&r5]1"),  # Cyclohexane ring
    )
    
    found_substructures = any(mol.HasSubstructMatch(pattern) for pattern in substructure_patterns)
    
    # Check for common monoterpenoid skeleton
    skeleton_pattern = Chem.MolFromSmarts("[C&r6,C&r5][C&r6,C&r5][C&r6,C&r5][C&r6,C&r5][C&r6,C&r5]1[C&r6,C&r5][C&r6,C&r5][C&r6,C&r5][C&r6,C&r5]1")
    has_skeleton = mol.HasSubstructMatch(skeleton_pattern)
    
    # Check for negative evidence (long hydrocarbon chains)
    negative_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    has_long_chain = mol.HasSubstructMatch(negative_pattern)
    
    if found_substructures and has_skeleton and not has_long_chain:
        return True, "Contains characteristic monoterpenoid substructures and skeleton"
    else:
        return False, "Does not exhibit typical monoterpenoid structural features"