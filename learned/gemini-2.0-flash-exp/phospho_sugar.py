"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is a monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General monosaccharide pattern: multiple hydroxyl groups and at least 3 carbons.
    # Count the number of carbons and oxygens with single bonds.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 2 )
    
    if c_count < 3 or o_count < 2: #Monosaccharide should have 3 carbons at least and 2 OH groups
        return False, "Not a monosaccharide: too few carbons or hydroxyl groups"
    
    # Phosphate group ester bond pattern (C-O-P=O)
    phosphate_ester_pattern = Chem.MolFromSmarts("[CX4,CX3][OX2][P](=[OX1])([OX2])([OX2])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_ester_pattern)

    # Cyclic phosphate group pattern (C-O-P-O-C)
    cyclic_phosphate_pattern = Chem.MolFromSmarts("[CX4,CX3][OX2][P]([OX2])([OX2])[OX2][CX4,CX3]")
    cyclic_phosphate_matches = mol.GetSubstructMatches(cyclic_phosphate_pattern)

    if not phosphate_matches and not cyclic_phosphate_matches:
        return False, "No phosphate group found directly linked via ester bond"
    
    #Check if there are enough O for the phospate group.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and any(x.GetAtomicNum() == 15 for x in atom.GetNeighbors()))
    if p_count > 0 and o_p_count < p_count * 3:
        return False, "Not enough oxygens around the phosphorus atoms"


    return True, "Contains a monosaccharide with a phosphate group attached via ester bond"