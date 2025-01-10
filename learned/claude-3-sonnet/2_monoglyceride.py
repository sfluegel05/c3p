"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone with a single fatty acid chain 
    attached via an ester bond at position 2, and free hydroxyls at positions 1 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens and kekulize
    mol = Chem.AddHs(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)

    # Core glycerol backbone pattern with ester at position 2
    # [OH]-CH2-CH(OR)-CH2-[OH]
    glycerol_pattern = Chem.MolFromSmarts("[OX2H]-[CH2X4]-[CHX4]-[CH2X4]-[OX2H]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with terminal hydroxyls found"
        
    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for 2-position attachment
    # Matches central carbon with ester attachment and terminal CH2OH groups
    position2_pattern = Chem.MolFromSmarts("[CH2X4][OX2H]-[CHX4]([OX2][CX3]=[OX1])-[CH2X4][OX2H]")
    if not mol.HasSubstructMatch(position2_pattern):
        return False, "Ester group not at position 2 or missing terminal hydroxyls"

    # Count carbons in the main chain after the ester
    # This pattern follows carbons after the ester carbonyl
    chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No fatty acid chain found"

    # Count primary alcohols (should be exactly 2)
    primary_oh_pattern = Chem.MolFromSmarts("[CH2X4][OX2H]")
    oh_matches = mol.GetSubstructMatches(primary_oh_pattern)
    if len(oh_matches) != 2:
        return False, f"Found {len(oh_matches)} primary hydroxyl groups, need exactly 2"

    # Basic size check
    if len(mol.GetAtoms()) < 10:  # Minimum size for smallest possible 2-monoglyceride
        return False, "Molecule too small to be a 2-monoglyceride"

    return True, "Contains glycerol backbone with single fatty acid chain at position 2"