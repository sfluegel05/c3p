"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine (CHEBI:89834)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    Must have a glycerol backbone with two acyl groups at positions 1 and 2, and a phosphoserine group at position 3.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule matches criteria, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for required components: glycerol backbone, two esters, phosphoserine group
    
    # 1. Find glycerol backbone with stereochemistry (sn-3 configuration)
    # Pattern matches [C@H] for chiral center with specific substitution pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[C@H](-[OX2])-[CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No proper glycerol backbone with sn-3 configuration"
    
    # 2. Check for two ester groups (acyl chains at positions 1 and 2)
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2]-C(=O)"))
    if len(ester_matches) < 2:
        return False, f"Found only {len(ester_matches)} ester groups, need at least 2"
    
    # 3. Verify phosphoserine group at position 3
    # Pattern matches phosphate connected to serine (NH2 and COOH groups)
    phosphoserine_pattern = Chem.MolFromSmarts("[OX2]-P(=O)([OX1-])-[OX2]-C(-[NH2])-C(=O)O")
    if not mol.HasSubstructMatch(phosphoserine_pattern):
        return False, "Missing phosphoserine group"
    
    # 4. Verify connectivity - phosphate must be attached to glycerol's 3rd carbon
    # Get atoms matching glycerol pattern
    matches = mol.GetSubstructMatches(glycerol_pattern)
    for match in matches:
        middle_carbon = match[1]  # [CH2]-[C@H]-[CH2] structure
        # Find attached oxygen atoms
        oxygen_neighbors = [n for n in mol.GetAtomWithIdx(middle_carbon).GetNeighbors()
                           if n.GetSymbol() == "O"]
        # One oxygen should be part of phosphate group
        phosphate_found = any(atom.GetSymbol() == "P" 
                            for o in oxygen_neighbors 
                            for atom in o.GetNeighbors() 
                            if atom.GetSymbol() == "P")
        if phosphate_found:
            return True, "Contains sn-glycerol backbone with two acyl groups and phosphoserine at position 3"
    
    return False, "Phosphate group not properly attached to glycerol backbone"