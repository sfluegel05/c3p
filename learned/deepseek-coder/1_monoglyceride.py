"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI:75541 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride has a glycerol backbone with a single fatty acid chain attached at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern with a single ester at position 1
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")

    # Check for the presence of a glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Find all ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Verify that the ester is attached to the first carbon of the glycerol backbone
    ester_carbon = ester_matches[0][0]  # Carbonyl carbon of the ester
    glycerol_atoms = mol.GetSubstructMatch(glycerol_pattern)
    first_carbon = glycerol_atoms[0]

    # Check if the ester carbon is bonded to the first carbon of the glycerol backbone
    ester_bonded_to_first_carbon = False
    for neighbor in mol.GetAtomWithIdx(ester_carbon).GetNeighbors():
        if neighbor.GetIdx() == first_carbon:
            ester_bonded_to_first_carbon = True
            break

    if not ester_bonded_to_first_carbon:
        return False, "Ester group not attached to the first carbon of the glycerol backbone"

    # Check for a fatty acid chain (long carbon chain attached to the ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds to verify the chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be a fatty acid"

    # Check molecular weight - 1-monoglycerides typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for 1-monoglyceride"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for 1-monoglyceride"
    if o_count != 4:
        return False, "Must have exactly 4 oxygens (1 ester group and 3 hydroxyls)"

    return True, "Contains glycerol backbone with a single fatty acid chain attached at position 1"