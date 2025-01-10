"""
Classifies: CHEBI:37739 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is any glycerolipid having a phosphate group ester-linked to a terminal carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol pattern (C-C-C with oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Detect phosphate group attached to glycerol
    phosphate_group_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_group_pattern)
    if not phosphate_matches:
        return False, "Missing phosphate group"

    # Check for ester-linked fatty acid chains
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;H2]-")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2 for fatty acid chains"

    # Check molecular weight to ensure it's in the typical range for glycerophospholipids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycerophospholipid"
    
    if mol_wt > 1000:
        return False, "Molecular weight too high for typical glycerophospholipid"

    return True, "Valid glycerophospholipid structure: glycerol backbone with attached phosphate and fatty acid chains"