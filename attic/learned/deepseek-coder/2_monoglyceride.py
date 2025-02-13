"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: CHEBI:75545 2-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride has a glycerol backbone with a fatty acid chain attached at position 2.

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

    # Look for glycerol backbone pattern (C-C-C with optional hydroxyl groups at positions 1 and 3)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester group found"

    # Check if the ester group is attached to the middle carbon of the glycerol backbone
    glycerol_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() >= 2]
    ester_attached_to_glycerol = False
    for match in ester_matches:
        ester_carbon = match[1]  # Carbon in the ester group
        if ester_carbon in glycerol_carbons:
            ester_attached_to_glycerol = True
            break
    if not ester_attached_to_glycerol:
        return False, "Ester group not attached to glycerol backbone"

    # Check for fatty acid chain (long carbon chain attached to ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Chain too short to be a fatty acid"

    # Check molecular weight - 2-monoglycerides typically >150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for 2-monoglyceride"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 6:
        return False, "Too few carbons for 2-monoglyceride"
    if o_count < 3:
        return False, "Must have at least 3 oxygens (1 ester group and 2 hydroxyl groups)"

    return True, "Contains glycerol backbone with a fatty acid chain attached at position 2 via an ester bond"