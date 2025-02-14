"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    A 1-O-acylglycerophosphoethanolamine has a glycerol backbone with a phosphate group
    attached to the 3-position, esterified with ethanolamine and a fatty acid at the 1-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for the glycerol backbone with substitution pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4](O)[CH2X4]OP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(glycerol_pattern):
          return False, "Glycerol backbone with phosphate and ethanolamine not found"
      
    # Check for ester group at position 1 of glycerol backbone
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH2X4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) != 1:
        return False, "Multiple/No glycerol substructure matches found"

    found_ester_1_position = False
    for match in ester_matches:
       ester_oxygen = match[0]
       for glycerol_match in glycerol_matches[0]:
        glycerol_mol = Chem.MolFromSmarts("[CH2X4][CHX4](O)[CH2X4]")
        glycerol_atoms = mol.GetSubstructMatches(glycerol_mol)[0]
        if mol.GetBondBetweenAtoms(glycerol_atoms[0],ester_oxygen) is not None:
            found_ester_1_position = True
            break
       if found_ester_1_position:
            break

    if not found_ester_1_position:
        return False, "No ester at position 1 of glycerol backbone"
   
    # Check for fatty acid chain (long carbon chain) at the 1-position
    # Looking for at least 4 carbons attached to the ester at position 1
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    found_fatty_chain = False
    for match in fatty_acid_matches:
         for glycerol_match in glycerol_matches[0]:
            glycerol_mol = Chem.MolFromSmarts("[CH2X4][CHX4](O)[CH2X4]")
            glycerol_atoms = mol.GetSubstructMatches(glycerol_mol)[0]
            if mol.GetBondBetweenAtoms(glycerol_atoms[0], mol.GetBondBetweenAtoms(match[0], match[1]).GetBeginAtomIdx()) is not None:
                found_fatty_chain = True
                break
         if found_fatty_chain:
            break

    if not found_fatty_chain:
         return False, "No fatty acid chain attached to ester at position 1"
    
    # Count rotatable bonds (for fatty acid chain check)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Fatty acid chain too short"

    # Check for at least one phosphorus and 3 oxygens (phosphate)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if p_count != 1:
        return False, "Must have exactly one phosphorus"
    if o_count < 6:
        return False, "Must have at least 6 oxygens, includes oxygens in glycerol, ester, phosphate and ethanolamine"
      
    return True, "Contains glycerol backbone with phosphate and ethanolamine and a fatty acid chain at the 1-position"