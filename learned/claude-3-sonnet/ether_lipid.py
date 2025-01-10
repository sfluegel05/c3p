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

    # Look for glycerol backbone patterns
    glycerol_patterns = [
        # Basic glycerol backbone
        Chem.MolFromSmarts("[OX2,OH1]-[CH2X4]-[CHX4]-[CH2X4]-[OX2,OH1]"),
        # Phospholipid glycerol backbone
        Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]-[CH2X4]-[CHX4]-[CH2X4]-[OX2]"),
        # More flexible glycerol pattern
        Chem.MolFromSmarts("[OX2,OH1]-[CH2X4]-[CH1,2X4]-[CH2X4]-[OX2,OH1]"),
        # Archaeal lipid glycerol pattern
        Chem.MolFromSmarts("[OX2]-[CH2X4]-[CHX4]-[CH2X4]-[OX2]"),
        # Phosphatidylcholine pattern
        Chem.MolFromSmarts("[CH3][N+]([CH3])([CH3])CCO[P](=[O])([O-])OC[CH2][CH]([CH2]O)")
    ]
    
    has_glycerol = False
    for pattern in glycerol_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_glycerol = True
            break
            
    if not has_glycerol:
        return False, "No glycerol backbone found"

    # Define ether patterns
    ether_patterns = [
        # Basic ether pattern
        Chem.MolFromSmarts("[CH2X4]-[OX2]-[CH2X4,CHX4]"),
        # Vinyl ether (plasmalogen)
        Chem.MolFromSmarts("[CH2X4]-[OX2]-[CH]=[CH]"),
        # Archaeal ether pattern
        Chem.MolFromSmarts("[CH2X4]-[OX2]-[CH2X4]-[CH2X4]"),
        # Branched ether pattern
        Chem.MolFromSmarts("[CH2X4]-[OX2]-[CH]([CH3])[CH2X4]")
    ]

    # Count ether linkages
    ether_count = 0
    for pattern in ether_patterns:
        if pattern is not None:
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                # Check if the ether is connected to glycerol backbone
                atom = mol.GetAtomWithIdx(match[0])
                if not atom.IsInRing():  # Exclude cyclic ethers
                    ether_count += 1

    if ether_count == 0:
        return False, "No ether linkages found"

    # Check for alkyl chains
    chain_patterns = [
        # Minimum 2-carbon chain
        Chem.MolFromSmarts("[CX4,CX3]-[CX4,CX3]"),
        # Branched chain
        Chem.MolFromSmarts("[CX4,CX3](-[CX4,CX3])-[CX4,CX3]")
    ]
    
    has_chain = False
    for pattern in chain_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_chain = True
            break

    if not has_chain:
        return False, "No alkyl chains found"

    # Basic composition checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 5:  # Minimum carbon requirement
        return False, "Too few carbons for ether lipid"
    if o_count < 2:  # Minimum oxygen requirement
        return False, "Too few oxygens for ether lipid"

    # Structural feature detection
    has_phosphate = mol.HasSubstructMatch(Chem.MolFromSmarts("[PX4](=O)([O-,OH])([O-,OH])"))
    has_ester = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4]-[OX2]-[CX3](=[OX1])[CX4]"))
    
    # Build classification reason
    features = []
    if has_phosphate:
        features.append("phosphate group")
    if has_ester:
        features.append("ester linkage")
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH2X4]-[OX2]-[CH]=[CH]")):
        features.append("vinyl ether")
    
    feature_str = ", ".join(features)
    if feature_str:
        feature_str = f" with {feature_str}"
    
    return True, f"Contains glycerol backbone with ether linkage to alkyl chain{feature_str}"