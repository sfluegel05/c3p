"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for glycerol backbone
    glycerol_smarts = "[C@@H](CO[P](=O)(O)O)O[C@@H](COC(=O)[#6])OC(=O)[#6]"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    if glycerol_pattern is None:
        return False, "Invalid glycerol SMARTS pattern"

    # Check for glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Define SMARTS pattern for inositol ring attached via phosphate to glycerol
    # Inositol ring: six-membered ring with six hydroxyl groups
    inositol_smarts = "C1(O)C(O)C(O)C(O)C(O)C1O"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if inositol_pattern is None:
        return False, "Invalid inositol SMARTS pattern"

    # Check for inositol ring
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Find the inositol ring atoms
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"

    # Identify phosphate groups attached to the inositol hydroxyls
    phosphate_count = 0
    for match in inositol_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # Check for oxygen atoms connected to the carbon atoms in the ring
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    for nei in neighbor.GetNeighbors():
                        if nei.GetAtomicNum() == 15:  # Phosphorus
                            phosphate_count += 1
                            break

    # Subtract the linkage phosphate to glycerol (already accounted in glycerol backbone)
    phosphate_count -= 1  # Adjust for the linkage phosphate

    if phosphate_count >= 1:
        return True, f"Phosphoinositide with {phosphate_count} phosphate group(s) on inositol ring"
    else:
        return False, "No phosphate groups on inositol ring beyond linkage phosphate"