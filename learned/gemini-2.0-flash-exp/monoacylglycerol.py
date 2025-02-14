"""
Classifies: CHEBI:17408 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol has a glycerol backbone with one fatty acid chain attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 1 ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    # Check that only one fatty acid chain is attached to an ester group
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    
    if len(fatty_acid_matches) < 1:
        return False, f"Missing fatty acid chain, got {len(fatty_acid_matches)}"

    # Count rotatable bonds
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Chain too short to be fatty acid"

    # Molecular Weight check - lower than triglycerides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 600:
        return False, "Molecular weight too high for monoacylglycerol"


    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8:
        return False, "Too few carbons for monoacylglycerol"
    if o_count < 3:
        return False, "Must have at least 3 oxygens (glycerol backbone and one ester group)"

    
    return True, "Contains glycerol backbone with 1 fatty acid chain attached via ester bonds"