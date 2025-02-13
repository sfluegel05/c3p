"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: CHEBI:17855 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D based on its SMILES string.
    Vitamin D compounds are characterized by a seco-steroid backbone with a conjugated triene system and hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the seco-steroid backbone with a conjugated triene system
    seco_steroid_pattern = Chem.MolFromSmarts("[C@H]1CC[C@@H]2[C@@]1(CCC/C2=C/C=C3/C[C@@H](O)CCC3=C)")
    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "No seco-steroid backbone with conjugated triene system found"

    # Check for at least one hydroxyl group (typical in vitamin D compounds)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl groups found"

    # Check for a flexible side chain at position 17
    side_chain_pattern = Chem.MolFromSmarts("[C@@H](CCCC(C)C)")
    if not mol.HasSubstructMatch(side_chain_pattern):
        return False, "No flexible side chain at position 17 found"

    # Check molecular weight - vitamin D compounds typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for vitamin D"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for vitamin D"
    if o_count < 1:
        return False, "Must have at least one oxygen (hydroxyl group)"

    return True, "Contains seco-steroid backbone with conjugated triene system and hydroxyl groups"