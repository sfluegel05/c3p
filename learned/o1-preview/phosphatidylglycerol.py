"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (C-C-C with hydroxyl groups)
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Check for two ester bonds attached to glycerol backbone
    ester_bond_pattern = Chem.MolFromSmarts("C(OC(=O)C)COC(=O)C")
    ester_bonds = mol.GetSubstructMatches(ester_bond_pattern)
    if not ester_bonds:
        return False, "Ester bonds to fatty acids not found on glycerol backbone"

    # Check for phosphate group attached to glycerol backbone
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group attached to glycerol backbone"

    # Check for additional glycerol attached via phosphodiester bond
    glycerol_phosphate_pattern = Chem.MolFromSmarts("COP(=O)(OCC(O)CO)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol connected via phosphodiester bond found"
    
    # Check for long fatty acid chains (at least 8 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)CCCCCCCC")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, "Less than two long fatty acid chains found"

    # Optionally, check molecular weight to ensure it's in the expected range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylglycerol"
    
    return True, "Molecule matches phosphatidylglycerol structure"