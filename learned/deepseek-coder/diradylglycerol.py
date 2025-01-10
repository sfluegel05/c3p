"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:17855 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol backbone with two substituent groups (acyl, alkyl, or alk-1-enyl) attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    
    # Get the glycerol backbone atoms
    glycerol_atoms = set(glycerol_matches[0])
    
    # Count substituent groups attached to the glycerol backbone
    substituent_count = 0
    for atom in mol.GetAtoms():
        if atom.GetIdx() in glycerol_atoms:
            continue  # Skip glycerol backbone atoms
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() in glycerol_atoms:
                # Check if the substituent is an acyl, alkyl, or alk-1-enyl group
                if atom.GetAtomicNum() == 6:  # Carbon
                    substituent_count += 1
                    break
    
    if substituent_count != 2:
        return False, f"Found {substituent_count} substituent groups attached to glycerol backbone, need exactly 2"

    # Check for long carbon chains attached to the substituents
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing long carbon chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Chains too short to be fatty acids or alkyl groups"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for diradylglycerol"
    if o_count < 3:
        return False, "Must have at least 3 oxygens (glycerol backbone and substituents)"

    return True, "Contains glycerol backbone with 2 substituent groups (acyl, alkyl, or alk-1-enyl)"