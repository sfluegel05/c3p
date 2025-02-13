"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:35781 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol is a glycerol backbone with three substituents - one acyl, one alkyl, and one alk-1-enyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C(O)-C(O)(O) with specific stereochemistry)
    glycerol_pattern = Chem.MolFromSmarts("[C@H]([O])([C@H]([O])[C@@H]([O]))")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found with correct stereochemistry"
    
    # Look for 3 substituents (acyl, alkyl, and alk-1-enyl)
    acyl_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    alkyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    alkenyl_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")
    
    substituents = []
    for match in matches[0]:
        atom = mol.GetAtomWithIdx(match)
        for neighbor in atom.GetNeighbors():
            if mol.HasSubstructMatch(acyl_pattern, atomIds=[neighbor.GetIdx()]):
                substituents.append("acyl")
            elif mol.HasSubstructMatch(alkyl_pattern, atomIds=[neighbor.GetIdx()]):
                substituents.append("alkyl")
            elif mol.HasSubstructMatch(alkenyl_pattern, atomIds=[neighbor.GetIdx()]):
                substituents.append("alk-1-enyl")
    
    # Check if there is exactly one acyl, one alkyl, and one alk-1-enyl substituent
    if len(set(substituents)) != 3 or len(substituents) != 3:
        return False, f"Found {len(substituents)} substituents, need exactly 3 (acyl, alkyl, alk-1-enyl)"
    
    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for triradylglycerol"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for triradylglycerol"
    
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 6:
        return False, "Must have exactly 6 oxygens"
    
    # Additional check for long carbon chains
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not any(mol.HasSubstructMatch(long_chain_pattern, atomIds=[neighbor.GetIdx()]) for match in matches[0] for neighbor in mol.GetAtomWithIdx(match).GetNeighbors()):
        return False, "No long carbon chains found"
    
    return True, "Glycerol backbone with one acyl, one alkyl, and one alk-1-enyl substituent"