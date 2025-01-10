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

    # Define SMARTS pattern for glycerophospholipid backbone
    glycerophospholipid_smarts = """
    [C;!R]-C(O[P](=O)(O)O)-O-C(=O)-[*]
    """
    glycerophospholipid_pattern = Chem.MolFromSmarts(glycerophospholipid_smarts)
    if glycerophospholipid_pattern is None:
        return False, "Invalid glycerophospholipid SMARTS pattern"

    # Check for glycerophospholipid backbone
    if not mol.HasSubstructMatch(glycerophospholipid_pattern):
        return False, "No glycerophospholipid backbone found"

    # Define SMARTS pattern for inositol ring attached via phosphate to glycerol
    # Inositol ring: six-membered ring with six hydroxyl groups
    inositol_smarts = "C1(O)C(O)C(O)C(O)C(O)C1O"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    if inositol_pattern is None:
        return False, "Invalid inositol SMARTS pattern"

    # Check for phosphatidylinositol structure (glycerol linked to inositol via phosphate)
    phosphatidylinositol_smarts = """
    [O;H0]-P(=O)([O;H0])[O]-[C@H]1([O])C(O)C(O)C(O)C(O)C1O
    """
    phosphatidylinositol_pattern = Chem.MolFromSmarts(phosphatidylinositol_smarts)
    if phosphatidylinositol_pattern is None:
        return False, "Invalid phosphatidylinositol SMARTS pattern"

    if not mol.HasSubstructMatch(phosphatidylinositol_pattern):
        return False, "No phosphatidylinositol structure found"

    # Find the inositol ring atoms
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"

    # Count the number of phosphate groups attached to inositol hydroxyls
    num_phosphates_on_inositol = 0

    for match in inositol_matches:
        inositol_atoms = set(match)
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon in inositol ring
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:  # Oxygen
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            # Check if oxygen is connected to a phosphate group
                            for nei in neighbor.GetNeighbors():
                                if nei.GetAtomicNum() == 15 and nei.GetIdx() not in inositol_atoms:
                                    num_phosphates_on_inositol += 1
                                    break

    # Exclude the linkage phosphate to glycerol backbone
    num_additional_phosphates = num_phosphates_on_inositol - 1  # Subtract linkage phosphate

    if num_additional_phosphates >= 1:
        return True, f"Phosphoinositide with {num_additional_phosphates} additional phosphate group(s) on inositol ring"
    else:
        return False, "No additional phosphate groups on inositol ring"