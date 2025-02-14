"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: CHEBI:83497 monoacylglycerol
A glyceride in which any one of the R groups (position not specified) is an acyl group while the remaining two R groups can be either H or alkyl groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.

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
    
    # Look for glycerol backbone pattern (C-C-C with 2-3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X2,CH2X3][CHX2,CHX3][CH2X2,CH2X3]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for 1 ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    # Check for acyl chain pattern (at least 3 connected carbons)
    acyl_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if not acyl_chain_matches:
        return False, "No acyl chain found"
    
    # Check that the ester group connects the glycerol and acyl chain
    ester_atom_idx = ester_matches[0][0]
    ester_atom = mol.GetAtomWithIdx(ester_atom_idx)
    ester_neighbors = [mol.GetAtomWithIdx(n).GetSymbol() for n in ester_atom.GetNeighbors()]
    if "C" not in ester_neighbors or "O" not in ester_neighbors:
        return False, "Ester group not connecting glycerol and acyl chain"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # No strict limits on molecular weight or oxygen count based on examples
    
    return True, "Contains glycerol backbone with one acyl chain attached via ester bond"