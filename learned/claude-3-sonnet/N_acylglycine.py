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

    # More flexible pattern for N-acylglycine core
    # Matches: R-C(=O)-N(R?)-CH2-C(=O)-O[H,Na,K,-]
    pattern = Chem.MolFromSmarts('[#6,#7,#8:1][C:2](=[O:3])[N:4]([#1,C:5])[CH2:6][C:7](=[O:8])[O:9][H,Na,K,-]')
    
    if not mol.HasSubstructMatch(pattern):
        return False, "No N-acylglycine moiety found"
    
    matches = mol.GetSubstructMatches(pattern)
    
    # Pattern to detect peptide bonds
    peptide_pattern = Chem.MolFromSmarts('[NH1,NH2][C:1](=[O:2])[C:3]([NH1,NH2])')
    
    # Pattern to detect carboxylic acid derivatives (esters, amides)
    deriv_pattern = Chem.MolFromSmarts('[C:1](=[O:2])[O,N][!H]')
    
    for match in matches:
        acyl_c = mol.GetAtomWithIdx(match[1])
        nitrogen = mol.GetAtomWithIdx(match[3])
        ch2 = mol.GetAtomWithIdx(match[5])
        acid_c = mol.GetAtomWithIdx(match[6])
        acid_o = mol.GetAtomWithIdx(match[8])
        
        # Check if this is part of a peptide chain
        if mol.HasSubstructMatch(peptide_pattern):
            peptide_matches = mol.GetSubstructMatches(peptide_pattern)
            is_peptide = False
            for pmatch in peptide_matches:
                if acid_c.GetIdx() == pmatch[0]:
                    is_peptide = True
                    break
            if is_peptide:
                continue
        
        # Verify CH2 group connectivity
        ch2_neighbors = list(ch2.GetNeighbors())
        if len(ch2_neighbors) != 2:
            continue
            
        # Check that CH2 neighbors are only N and acid C
        ch2_neighbor_ids = {n.GetIdx() for n in ch2_neighbors}
        if not ch2_neighbor_ids == {nitrogen.GetIdx(), acid_c.GetIdx()}:
            continue
        
        # Check carboxylic acid group isn't derivatized
        acid_deriv_matches = mol.GetSubstructMatches(deriv_pattern)
        is_derivatized = False
        for dmatch in acid_deriv_matches:
            if acid_c.GetIdx() == dmatch[0]:
                # Check if it's not just a salt or charged form
                if not any(a.GetFormalCharge() < 0 for a in acid_c.GetNeighbors()):
                    is_derivatized = True
                    break
        if is_derivatized:
            continue
            
        # Accept both unsubstituted and N-substituted glycines
        if nitrogen.GetDegree() in [2, 3]:  # NH or N-substituted
            return True, "Contains N-acylglycine structure: R-C(=O)-N(R?)-CH2-C(=O)-OH"
    
    return False, "Does not contain valid N-acylglycine structure"