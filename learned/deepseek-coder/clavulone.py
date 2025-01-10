"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: Clavulone (a class of esterified prostanoids obtained from marine corals)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Clavulones are esterified prostanoids with a cyclopentane ring, hydroxyl group, and ester groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cyclopentane ring with a ketone group (C=O)
    cyclopentane_ketone_pattern = Chem.MolFromSmarts("[C]1[C][C][C][C]1(=O)")
    if not mol.HasSubstructMatch(cyclopentane_ketone_pattern):
        return False, "No cyclopentane ring with ketone group found"

    # Check for hydroxyl group (OH) attached to the cyclopentane ring
    hydroxyl_pattern = Chem.MolFromSmarts("[C][C]([OH])[C]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group attached to cyclopentane ring found"

    # Check for ester groups (acetates or similar)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found"

    # Check for long hydrocarbon chains with double bonds (C=C)
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    hydrocarbon_chain_matches = mol.GetSubstructMatches(hydrocarbon_chain_pattern)
    if len(hydrocarbon_chain_matches) < 1:
        return False, "No long hydrocarbon chains with double bonds found"

    # Check for halogens (Cl, Br, I) if present
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I]")
    halogen_matches = mol.GetSubstructMatches(halogen_pattern)
    if len(halogen_matches) > 0:
        return True, "Contains cyclopentane ring with ketone, hydroxyl, ester groups, and halogen"
    else:
        return True, "Contains cyclopentane ring with ketone, hydroxyl, and ester groups"