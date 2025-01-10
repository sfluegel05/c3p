"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol is a glycerophosphoinositol having one phosphatidyl group
    esterified to one of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for diacylglycerol backbone with esterified fatty acids at sn-1 and sn-2 positions
    diacylglycerol_smarts = """
    [C@@H]([CH2]OC(=O)C)COC(=O)C
    """
    diacylglycerol_pattern = Chem.MolFromSmarts(diacylglycerol_smarts)
    if not mol.HasSubstructMatch(diacylglycerol_pattern):
        return False, "No diacylglycerol backbone with esterified fatty acids found"

    # Check for glycerol phosphate group attached at sn-3 position
    glycerol_phosphate_smarts = """
    [C@@H](CO[P](=O)(O)O)([CH2]OC(=O)C)COC(=O)C
    """
    glycerol_phosphate_pattern = Chem.MolFromSmarts(glycerol_phosphate_smarts)
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No phosphate group attached to glycerol backbone found"

    # Check for inositol ring attached via phosphodiester bond
    inositol_ring_smarts = """
    O[P](=O)(OC1C(O)C(O)C(O)C(O)C1O)O
    """
    inositol_ring_pattern = Chem.MolFromSmarts(inositol_ring_smarts)
    if not mol.HasSubstructMatch(inositol_ring_pattern):
        return False, "No inositol ring attached via phosphodiester bond found"

    # Ensure the phosphate group connects glycerol and inositol correctly
    # Find phosphate atoms
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    phosphate_found = False
    for match in phosphate_matches:
        phosphorous_atom_idx = match[1]  # Index of the phosphorous atom
        phosphorous_atom = mol.GetAtomWithIdx(phosphorous_atom_idx)
        
        # Get neighbors of the phosphorous atom
        neighbors = phosphorous_atom.GetNeighbors()
        
        # Check if one neighbor is connected to glycerol and another to inositol
        connected_to_glycerol = False
        connected_to_inositol = False
        
        for neighbor in neighbors:
            neighbor_idx = neighbor.GetIdx()
            atom_env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, neighbor_idx)
            submol = Chem.PathToSubmol(mol, atom_env)
            
            if submol.HasSubstructMatch(diacylglycerol_pattern):
                connected_to_glycerol = True
            if submol.HasSubstructMatch(inositol_ring_pattern):
                connected_to_inositol = True

        if connected_to_glycerol and connected_to_inositol:
            phosphate_found = True
            break

    if not phosphate_found:
        return False, "Phosphate group not properly connecting glycerol and inositol"

    return True, "Molecule is a phosphatidylinositol"