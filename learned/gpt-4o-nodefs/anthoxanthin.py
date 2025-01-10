"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid, typically containing two aromatic rings
    connected by a three-carbon bridge forming a pyran ring, often with several hydroxyl
    groups and sometimes glycosidic attachments.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for flavonoid backbone (C6-C3-C6)
    flavonoid_pattern = Chem.MolFromSmarts("c1ccccc1C2=COC(=O)c3ccccc23")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "Does not contain typical flavonoid backbone structure"

    # Check for hydroxyl and ketone groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and any(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in atom.GetBonds()))
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if hydroxyl_count < 3 or len(ketone_matches) < 1:
        return False, "Does not have sufficient hydroxyl or ketone groups"

    # Optional: Check for glycoside groups, identified by oxygen-containing rings
    glycoside_pattern = Chem.MolFromSmarts("OC1COC(O)C(O)C1O")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycoside groups, which are common but not essential"

    return True, "Contains anthoxanthin characteristic flavonoid backbone with appropriate functional groups"