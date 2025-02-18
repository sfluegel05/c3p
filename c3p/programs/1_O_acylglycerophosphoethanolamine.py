"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    These compounds have:
    - A glycerol backbone
    - An acyl group (ester) at sn-1 position
    - A hydroxyl at sn-2 position
    - A phosphoethanolamine group at sn-3 position
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (3-carbon chain with oxygen substituents)
    glycerol = Chem.MolFromSmarts("[OX2][CH2][CH]([OX2])[CH2][OX2]")
    if not mol.HasSubstructMatch(glycerol):
        return False, "Missing glycerol backbone"

    # Check for ester group
    ester = Chem.MolFromSmarts("[CX3](=O)[OX2][CH2]")
    if not mol.HasSubstructMatch(ester):
        return False, "Missing ester group"

    # Check for phosphate group (both neutral and charged forms)
    phosphate = Chem.MolFromSmarts("[OX2][P](=[O])([OX2])[OX2]")
    if not mol.HasSubstructMatch(phosphate):
        return False, "Missing phosphate group"

    # Check for ethanolamine group (both neutral and charged forms)
    ethanolamine = Chem.MolFromSmarts("[OX2]CC[NX3,NX4+]")
    if not mol.HasSubstructMatch(ethanolamine):
        return False, "Missing ethanolamine group"

    # Verify it's not a phosphocholine
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[NX4+](C)(C)(C)")):
        return False, "Contains quaternary amine (phosphocholine)"

    # Count key atoms to ensure correct composition
    p_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#15]")))
    if p_count != 1:
        return False, f"Should have exactly 1 phosphorus, found {p_count}"

    n_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#7]")))
    if n_count != 1:
        return False, f"Should have exactly 1 nitrogen, found {n_count}"

    # Count ester groups - should be exactly 1
    ester_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OX2]")))
    if ester_count != 1:
        return False, f"Should have exactly 1 ester group, found {ester_count}"

    # Check for minimum chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 7:  # Minimum: 3 (glycerol) + 2 (ethanolamine) + 2 (acyl)
        return False, "Carbon chain too short"

    # Verify specific connectivity: ester at sn-1, OH at sn-2, phosphate at sn-3
    complete_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CH2][CH]([OX2H,OX2-])[CH2][OX2]P")
    if not mol.HasSubstructMatch(complete_pattern):
        return False, "Incorrect arrangement of substituents on glycerol backbone"

    return True, "Contains glycerol backbone with phosphoethanolamine group and single acyl substituent at sn-1"