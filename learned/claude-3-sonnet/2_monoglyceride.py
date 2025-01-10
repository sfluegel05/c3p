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

    # Look for glycerol backbone pattern (C-C-C with oxygens)
    # More flexible pattern that doesn't require explicit hydrogens
    glycerol_pattern = Chem.MolFromSmarts("[OX2]-[CH2X4]-[CHX4]-[CH2X4]-[OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])-[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for 2-position attachment
    # More flexible pattern for central carbon with ester
    position2_pattern = Chem.MolFromSmarts("[OX2]-[CH2X4]-[CHX4](-[OX2][CX3]=[OX1])-[CH2X4]-[OX2]")
    if not mol.HasSubstructMatch(position2_pattern):
        return False, "Ester group not at position 2"

    # Count oxygens (should be 4: 2 from hydroxyls, 2 from ester)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 4:
        return False, f"Found {o_count} oxygens, need exactly 4"

    # Verify fatty acid chain length (at least 2 carbons after ester)
    chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[#6]-[#6]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Fatty acid chain too short"

    # Count total carbons (should be at least 6: 3 from glycerol + at least 3 from fatty acid)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Molecule too small to be a 2-monoglyceride"

    return True, "Contains glycerol backbone with single fatty acid chain at position 2"