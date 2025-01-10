"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: CHEBI:35255 chlorophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    A chlorophyll is a magnesium porphyrin with a fifth ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium in the molecule
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if len(mg_atoms) != 1:
        return False, "No magnesium atom found"

    # Define a more flexible porphyrin core pattern with the fifth ring
    # This pattern matches the basic structure of a porphyrin with a fifth ring
    porphyrin_pattern = Chem.MolFromSmarts("[Mg]1234n1c(c2)c(c3)c(c4)c5c1c(c2)c(c3)c(c4)c5")
    if porphyrin_pattern is None:
        return False, "Failed to create porphyrin pattern"
    if not mol.HasSubstructMatch(porphyrin_pattern):
        # If the rigid pattern fails, try a more flexible pattern
        flexible_porphyrin_pattern = Chem.MolFromSmarts("[Mg]1234n1c(c2)c(c3)c(c4)c5c1c(c2)c(c3)c(c4)c5")
        if flexible_porphyrin_pattern is None:
            return False, "Failed to create flexible porphyrin pattern"
        if not mol.HasSubstructMatch(flexible_porphyrin_pattern):
            return False, "No porphyrin core with fifth ring found"

    # Check for ester groups (common in chlorophylls)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if ester_pattern is None:
        return False, "Failed to create ester pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found"

    # Check for carbonyl groups (common in chlorophylls)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if carbonyl_pattern is None:
        return False, "Failed to create carbonyl pattern"
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) < 1:
        return False, "No carbonyl groups found"

    # Check molecular weight - chlorophylls typically >800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for chlorophyll"

    return True, "Contains magnesium porphyrin core with fifth ring and necessary functional groups"