"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid has a glycerol backbone with at least one alkyl chain 
    attached via an ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ether linkage (R-O-R') not part of an ester group
    # We'll look for oxygen connected to two carbons (ether)
    ether_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    if not ether_matches:
        return False, "No ether linkages found"

    # Look for long alkyl chains (at least 4 carbons) connected to ether oxygen
    long_chain_pattern = Chem.MolFromSmarts("[CX4]-[CX4]-[CX4]-[CX4]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long alkyl chains found"

    # Count carbons and oxygens to ensure reasonable size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for ether lipid"
    if o_count < 3:
        return False, "Too few oxygens for ether lipid"

    # Additional checks for common ether lipid features
    
    # Check for phosphate group (common in phosphoether lipids)
    has_phosphate = mol.HasSubstructMatch(Chem.MolFromSmarts("[PX4](=O)([O-,OH])([O-,OH])"))
    
    # Check for ester groups
    has_ester = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4]-[OX2]-[CX3](=[OX1])[CX4]"))
    
    # Check for vinyl ether (plasmalogen) pattern
    has_vinyl_ether = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4]-[OX2]-[CX3]=[CX3]"))

    # Build detailed reason string
    features = []
    if has_phosphate:
        features.append("phosphate group")
    if has_ester:
        features.append("ester linkage")
    if has_vinyl_ether:
        features.append("vinyl ether")
    
    feature_str = ", ".join(features)
    if feature_str:
        feature_str = f" with {feature_str}"
    
    return True, f"Contains glycerol backbone with ether linkage to alkyl chain{feature_str}"