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

    # Look for glycerol backbone connected to proper substituents
    glycerol_patterns = [
        # Basic glycerol backbone
        Chem.MolFromSmarts("[OX2,OH1]-[CH2X4]-[CHX4]-[CH2X4]-[OX2,OH1]"),
        # Phospholipid glycerol backbone
        Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]-[CH2X4]-[CHX4]-[CH2X4]-[OX2]"),
        # More flexible glycerol pattern
        Chem.MolFromSmarts("[OX2,OH1]-[CH2X4]-[CH1,2X4]-[CH2X4]-[OX2,OH1]")
    ]
    
    has_glycerol = any(mol.HasSubstructMatch(pattern) for pattern in glycerol_patterns)
    if not has_glycerol:
        return False, "No glycerol backbone found"

    # Look for ether linkages specifically connected to glycerol backbone
    ether_patterns = [
        # Saturated ether on glycerol
        Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2]-[CH2X4,CHX3]-[!O;!N;!S]"),
        # Vinyl ether (plasmalogen)
        Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2]-[CH1X3]=[CHX3]"),
        # General ether pattern excluding esters
        Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2]-[CX4,CX3;!$(C=O)]")
    ]
    
    # Count ether linkages that are part of the glycerol moiety
    glycerol_ether_count = 0
    for pattern in ether_patterns:
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # Check if the ether is connected to the glycerol backbone
            if any(mol.GetAtomWithIdx(match[0]).IsInRing() for match in matches):
                continue
            glycerol_ether_count += 1

    if glycerol_ether_count == 0:
        return False, "No ether linkages found connected to glycerol backbone"

    # Look for alkyl chains (more flexible length requirement)
    chain_patterns = [
        # Minimum 4-carbon chain
        Chem.MolFromSmarts("[CX4,CX3]-[CX4,CX3]-[CX4,CX3]-[CX4,CX3]"),
        # Unsaturated chain
        Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    ]
    
    has_chain = any(mol.HasSubstructMatch(pattern) for pattern in chain_patterns)
    if not has_chain:
        return False, "No suitable alkyl chains found"

    # Basic size and composition checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8:  # Relaxed carbon requirement
        return False, "Too few carbons for ether lipid"
    if o_count < 3:
        return False, "Too few oxygens for ether lipid"

    # Check for specific structural features
    has_phosphate = mol.HasSubstructMatch(Chem.MolFromSmarts("[PX4](=O)([O-,OH])([O-,OH])"))
    has_ester = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4]-[OX2]-[CX3](=[OX1])[CX4]"))
    has_plasmalogen = mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2X4]-[OX2]-[CH1X3]=[CHX3]"))
    
    # Exclude molecules with too many rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 2:
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