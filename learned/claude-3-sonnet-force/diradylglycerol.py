"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:36764 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol backbone with two substituent groups (acyl, alkyl, or alk-1-enyl) attached at any two positions.

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

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X3]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester/ether groups (-O-C or -O-C=O)
    ester_ether_pattern = Chem.MolFromSmarts("[OX2][CX3,CX3]=O")
    ester_ether_matches = mol.GetSubstructMatches(ester_ether_pattern)
    if len(ester_ether_matches) != 2:
        return False, f"Found {len(ester_ether_matches)} ester/ether groups, need exactly 2"

    # Check for alkyl/acyl/alk-1-enyl chains (carbon chains attached to esters/ethers)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 2:
        return False, f"Missing substituent chains, got {len(chain_matches)}"

    # Count rotatable bonds to verify chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Substituent chains too short"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for diradylglycerol"
    if o_count != 4:
        return False, "Must have exactly 4 oxygens (2 ester/ether groups, 1 backbone)"

    return True, "Contains glycerol backbone with 2 substituent groups attached"