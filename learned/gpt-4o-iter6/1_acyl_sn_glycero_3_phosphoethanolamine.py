"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    This class is characterized as a 1-O-acylglycerophosphoethanolamine with (R)-configuration at sn-2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycerol backbone with chiral center at sn-2
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](O)C(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with (R)-configuration at sn-2 found"

    # Check for phosphate group linked in the PHOS motif
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphoethanolamine group found"
    
    # Check for ester linkage at sn-1
    ester_pattern = Chem.MolFromSmarts("O[C@H](O)COC(=O)")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage at sn-1 position found"

    # Ensure acyl chain presence
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)[C,C]")
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No acyl chain detected connected via ester"

    return True, "Molecule matches 1-acyl-sn-glycero-3-phosphoethanolamine structural criteria"

__metadata__ = { 
    'chemical_class': {   
        'id': None,
        'name': '1-acyl-sn-glycero-3-phosphoethanolamine',
        'definition': 'A 1-O-acylglycerophosphoethanolamine having (R)-configuration.'
    }
}