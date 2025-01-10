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

    # Pattern for alpha-amino acid core with any substituent on nitrogen
    # [C@H] ensures S stereochemistry (L-amino acid)
    amino_acid_pattern = Chem.MolFromSmarts("[C@H](N)(C(=O)[OH1,O-])[#6,#1]")
    if not mol.HasSubstructMatch(amino_acid_pattern, useChirality=True):
        return False, "No L-alpha-amino acid core found"

    # Pattern for N-acyl group (R-C=O-N)
    # The nitrogen should be connected to the alpha carbon
    n_acyl_pattern = Chem.MolFromSmarts("[#6,#1]C(=O)[NH1]C([#6,#1])(C(=O)[OH1,O-])[#6,#1]")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl group found"

    # Get all matches of the amino acid pattern
    matches = mol.GetSubstructMatches(amino_acid_pattern, useChirality=True)
    
    for match in matches:
        alpha_carbon_idx = match[0]
        nitrogen_idx = match[1]
        
        # Get the alpha carbon and check its chirality explicitly
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
        if alpha_carbon.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
            
        # Get the nitrogen atom and check its substitution
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Nitrogen should have exactly one hydrogen (NH) since it's N-acylated
        if nitrogen.GetTotalNumHs() != 1:
            continue
            
        # Check that nitrogen is connected to a carbonyl carbon
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                for n2 in neighbor.GetNeighbors():
                    if n2.GetAtomicNum() == 8 and n2.GetIsAromatic() == False:  # Oxygen
                        if sum(1 for b in neighbor.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE) == 1:
                            return True, "Contains N-acyl-L-alpha-amino acid structure"
    
    return False, "No valid N-acyl-L-alpha-amino acid structure found with correct substitution pattern"