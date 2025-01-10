"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:58080 alpha-amino acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid zwitterion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an alpha-amino acid zwitterion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for NH3+ group
    nh3_pattern = Chem.MolFromSmarts("[NH3+]")
    nh3_matches = mol.GetSubstructMatches(nh3_pattern)
    if not nh3_matches:
        return False, "No NH3+ group found"
    
    # Look for carboxylate group (C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[C](=[O])[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not carboxylate_matches:
        return False, "No carboxylate group found"
    
    # Pattern for alpha-amino acid zwitterion core structure:
    # [NH3+] connected to carbon which is connected to C(=O)[O-]
    core_pattern = Chem.MolFromSmarts("[NH3+][C;X4][C](=[O])[O-]")
    core_matches = mol.GetSubstructMatches(core_pattern)
    
    if not core_matches:
        return False, "No alpha-amino acid zwitterion core structure found"
    
    # For each potential core, verify it's a proper alpha-amino acid structure
    for match in core_matches:
        nh3_idx, alpha_c_idx, carboxyl_c_idx = match[0:3]
        
        # Get the alpha carbon atom
        alpha_carbon = mol.GetAtomWithIdx(alpha_c_idx)
        
        # Check if alpha carbon has at least one hydrogen
        # (total valence - explicit valence = number of implicit hydrogens)
        if alpha_carbon.GetTotalNumHs() > 0:
            return True, "Contains alpha-amino acid zwitterion structure with NH3+ and COO- groups"
    
    return False, "No valid alpha-amino acid zwitterion structure found"