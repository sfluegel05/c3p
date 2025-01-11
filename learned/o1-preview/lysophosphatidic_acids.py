"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    Lysophosphatidic acids are monoacylglycerol phosphates obtained by hydrolytic removal of one of the two acyl groups of any phosphatidic acid or derivatives therein.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    # Glycerol backbone with phosphate group attached at sn-3 position
    glycerol_phosphate_pattern = Chem.MolFromSmarts("OCC(O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol phosphate backbone found"
    
    # Identify glycerol backbone atoms
    glycerol_match = mol.GetSubstructMatch(glycerol_phosphate_pattern)
    glycerol_carbons = [atom_idx for atom_idx in glycerol_match if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6]
    
    # Define ester bond pattern (excluding phosphoester bonds)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Count ester bonds connected to glycerol carbons
    acyl_ester_count = 0
    for match in ester_matches:
        ester_o_idx = match[2]  # index of the oxygen in ester bond
        ester_c_idx = match[3]  # index of the carbon connected to ester oxygen (should be glycerol carbon)
        if ester_c_idx in glycerol_carbons:
            acyl_ester_count += 1
    
    if acyl_ester_count != 1:
        return False, f"Expected 1 acyl chain attached via ester bond to glycerol backbone, found {acyl_ester_count}"
    
    # Ensure there are no additional acyl chains
    # Count total ester bonds excluding phosphoester bonds
    total_acyl_esters = len(ester_matches)
    if total_acyl_esters > 1:
        return False, f"Found {total_acyl_esters} acyl ester bonds, expected only 1"
    
    return True, "Molecule is a lysophosphatidic acid with glycerol phosphate backbone and one acyl chain"