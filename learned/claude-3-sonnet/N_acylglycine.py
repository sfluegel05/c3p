"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N-acylglycine
An N-acyl-amino acid in which amino acid specified is glycine.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    N-acylglycines have a glycine moiety (NH-CH2-COOH) with an acyl group (R-C=O-) 
    attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N-acylglycine pattern:
    # [C:1](=[O:2])[N:3][CH2:4][C:5](=[O:6])[O,OH]
    # This matches:
    # - Any carbon with double-bonded oxygen (the acyl group)
    # - Connected to a nitrogen
    # - Connected to CH2
    # - Connected to a carboxyl group (handles both acid and charged forms)
    pattern = Chem.MolFromSmarts('[C:1](=[O:2])[N:3][CH2:4][C:5](=[O:6])[O,OH]')
    
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No N-acylglycine moiety found"
    
    # For each match, verify it's a true N-acylglycine:
    for match in matches:
        acyl_c = mol.GetAtomWithIdx(match[0])
        nitrogen = mol.GetAtomWithIdx(match[2])
        ch2 = mol.GetAtomWithIdx(match[3])
        acid_c = mol.GetAtomWithIdx(match[4])
        
        # Verify acyl carbon has a non-oxygen neighbor (R group)
        has_r_group = False
        for neighbor in acyl_c.GetNeighbors():
            if neighbor.GetAtomicNum() != 8 and neighbor.GetAtomicNum() != 7:
                has_r_group = True
                break
                
        if not has_r_group:
            continue
            
        # Verify CH2 is only connected to N and acid C
        if len(list(ch2.GetNeighbors())) != 2:
            continue
            
        # All checks passed
        return True, "Contains N-acylglycine structure: R-C(=O)-NH-CH2-C(=O)-O"
        
    return False, "Contains similar groups but not in N-acylglycine arrangement"