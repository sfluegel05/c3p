"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is a member of the group of substituted penams containing:
    - Two methyl substituents at position 2
    - A carboxylate substituent at position 3
    - A carboxamido group at position 6

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the penicillin core structure pattern
    penicillin_core = Chem.MolFromSmarts("[C]1([C@@H]2N([C@@H]1=O)[C@@H](C(S2)(C)C)C(=O)O)C")
    if not mol.HasSubstructMatch(penicillin_core):
        return False, "No penicillin core structure found"

    # Check for two methyl groups at position 2
    methyl_pattern = Chem.MolFromSmarts("[C]1([C@@H]2N([C@@H]1=O)[C@@H](C(S2)([CH3])([CH3]))C(=O)O)C")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) == 0:
        return False, "No methyl groups found at position 2"

    # Check for carboxylate group at position 3
    carboxylate_pattern = Chem.MolFromSmarts("[C]1([C@@H]2N([C@@H]1=O)[C@@H](C(S2)(C)C)C(=O)[O-])C")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) == 0:
        return False, "No carboxylate group found at position 3"

    # Check for carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts("[C]1([C@@H]2N([C@@H]1=O)[C@@H](C(S2)(C)C)C(=O)N)C")
    carboxamido_matches = mol.GetSubstructMatches(carboxamido_pattern)
    if len(carboxamido_matches) == 0:
        return False, "No carboxamido group found at position 6"

    # Additional check for molecular weight (penicillins are typically > 300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for penicillin"

    return True, "Contains penicillin core structure with required substituents"