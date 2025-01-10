"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid/carboxylate group
    carboxyl_pattern = Chem.MolFromSmarts("[$([CX3](=[OX1])[OX2H1]),$([CX3](=[OX1])[OX1-])]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid/carboxylate group found"

    # Look for alpha-amino acid core with N-acyl
    # More general pattern that matches both protonated and deprotonated forms
    # and different representations of the acyl group
    amino_acid_pattern = Chem.MolFromSmarts("""
        [C:1]([*,H])([N:2]([H,#6,#7,#8])[C:3](=[O:4])[#6,#7,#8,#16])[C:5](=[O:6])[O,O-]
    """)
    
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not matches:
        return False, "No N-acyl alpha-amino acid core structure found"

    # Check each match for correct configuration
    for match in matches:
        alpha_carbon_idx = match[0]
        nitrogen_idx = match[1]
        carbonyl_idx = match[2]
        
        # Get the alpha carbon and check its chirality
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
        
        # Skip if not chiral
        if alpha_carbon.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
            
        # Check for L configuration (S in CIP)
        # Need to check both CCW and CW because SMILES representation can vary
        chiral_tag = alpha_carbon.GetChiralTag()
        if chiral_tag in [Chem.ChiralType.CHI_TETRAHEDRAL_CCW, Chem.ChiralType.CHI_TETRAHEDRAL_CW]:
            # Verify N-acyl group
            nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
            carbonyl = mol.GetAtomWithIdx(carbonyl_idx)
            
            # Check nitrogen substitution
            if nitrogen.GetTotalNumHs() > 1:  # Should not be NH2
                continue
                
            # Check carbonyl is part of acyl group
            if carbonyl.GetTotalNumHs() > 0:  # Carbonyl carbon should have no H
                continue
                
            return True, "Contains N-acyl-L-alpha-amino acid structure"
            
    return False, "No valid N-acyl-L-alpha-amino acid structure found with correct stereochemistry"