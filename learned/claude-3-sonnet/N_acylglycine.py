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

    # Look for basic N-acylglycine pattern
    pattern = Chem.MolFromSmarts('[C:1](=[O:2])[N:3][CH2:4][C:5](=[O:6])[O,OH]')
    
    if not mol.HasSubstructMatch(pattern):
        return False, "No N-acylglycine moiety found"
    
    matches = mol.GetSubstructMatches(pattern)
    
    # Pattern to detect peptide bonds
    peptide_pattern = Chem.MolFromSmarts('[NX3:1][C:2](=[O:3])[C:4]([NX3:5])')
    
    for match in matches:
        acyl_c = mol.GetAtomWithIdx(match[0])
        nitrogen = mol.GetAtomWithIdx(match[2])
        ch2 = mol.GetAtomWithIdx(match[3])
        acid_c = mol.GetAtomWithIdx(match[4])
        
        # Check if nitrogen is part of a peptide chain
        # Get all atoms in peptide bonds
        peptide_matches = mol.GetSubstructMatches(peptide_pattern)
        peptide_nitrogens = set()
        for pmatch in peptide_matches:
            peptide_nitrogens.add(pmatch[0])
            peptide_nitrogens.add(pmatch[4])
        
        # Skip if this nitrogen is part of a peptide chain
        if nitrogen.GetIdx() in peptide_nitrogens:
            continue
            
        # Verify CH2 is only connected to N and acid C
        ch2_neighbors = list(ch2.GetNeighbors())
        if len(ch2_neighbors) != 2:
            continue
            
        # Check that CH2 neighbors are only N and acid C
        ch2_neighbor_ids = {n.GetIdx() for n in ch2_neighbors}
        if not ch2_neighbor_ids == {nitrogen.GetIdx(), acid_c.GetIdx()}:
            continue
        
        # Verify carboxylic acid group
        acid_neighbors = list(acid_c.GetNeighbors())
        has_oh = False
        for neighbor in acid_neighbors:
            if neighbor.GetAtomicNum() == 8:
                has_oh = True
                break
        if not has_oh:
            continue
            
        # If we get here, we have a valid N-acylglycine
        return True, "Contains N-acylglycine structure: R-C(=O)-NH-CH2-C(=O)-OH"
    
    return False, "Contains similar groups but not in correct N-acylglycine arrangement"