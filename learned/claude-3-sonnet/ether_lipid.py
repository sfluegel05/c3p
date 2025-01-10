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
    attached via an ether linkage (including vinyl ethers).

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

    # Look for glycerol backbone with specific substitution pattern
    # Match both regular and phospholipid glycerol backbones
    glycerol_patterns = [
        # Regular glycerol backbone with ether/ester positions
        Chem.MolFromSmarts("[OX2]-[CH2X4]-[CHX4]-[CH2X4]-[OX2]"),
        # Phospholipid glycerol backbone
        Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]-[CH2X4]-[CHX4]-[CH2X4]-[OX2]")
    ]
    
    has_glycerol = any(mol.HasSubstructMatch(pattern) for pattern in glycerol_patterns)
    if not has_glycerol:
        return False, "No glycerol backbone found"

    # Look for ether linkages (both saturated and vinyl)
    ether_patterns = [
        # Saturated ether
        Chem.MolFromSmarts("[CH2X4]-[OX2]-[CH2X4,CHX3]"),
        # Vinyl ether (plasmalogen)
        Chem.MolFromSmarts("[CH2X4]-[OX2]-[CH1X3]=[CHX3]"),
        # Other unsaturated ether variations
        Chem.MolFromSmarts("[CH2X4]-[OX2]-[CX3]")
    ]
    
    ether_matches = []
    for pattern in ether_patterns:
        matches = mol.GetSubstructMatches(pattern)
        ether_matches.extend(matches)
    
    if not ether_matches:
        return False, "No ether linkages found"

    # Look for long alkyl chains (at least 8 carbons) connected to ether oxygen
    long_chain_patterns = [
        # Saturated chain
        Chem.MolFromSmarts("[CX4]-[CX4]-[CX4]-[CX4]-[CX4]-[CX4]-[CX4]-[CX4]"),
        # Chain with unsaturation
        Chem.MolFromSmarts("[$([CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3])]")
    ]
    
    has_long_chain = any(mol.HasSubstructMatch(pattern) for pattern in long_chain_patterns)
    if not has_long_chain:
        return False, "No long alkyl chains found"

    # Verify the molecule has appropriate size and composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 12:  # Increased minimum carbon requirement
        return False, "Too few carbons for ether lipid"
    if o_count < 3:
        return False, "Too few oxygens for ether lipid"

    # Check for specific structural features
    has_phosphate = mol.HasSubstructMatch(Chem.MolFromSmarts("[PX4](=O)([O-,OH])([O-,OH])"))
    has_ester = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4]-[OX2]-[CX3](=[OX1])[CX4]"))
    has_plasmalogen = mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2X4]-[OX2]-[CH1X3]=[CHX3]"))
    
    # Exclude molecules with too many rings (likely not lipids)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 2:  # Allow up to 2 rings (some lipids have cyclic groups)
        return False, "Too many rings for ether lipid"

    # Build detailed reason string
    features = []
    if has_phosphate:
        features.append("phosphate group")
    if has_ester:
        features.append("ester linkage")
    if has_plasmalogen:
        features.append("vinyl ether (plasmalogen)")
    
    feature_str = ", ".join(features)
    if feature_str:
        feature_str = f" with {feature_str}"
    
    return True, f"Contains glycerol backbone with ether linkage to alkyl chain{feature_str}"