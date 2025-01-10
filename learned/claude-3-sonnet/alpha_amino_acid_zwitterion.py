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

    # Check for counter ions like Na+, K+, etc.
    metal_pattern = Chem.MolFromSmarts("[Na+,K+,Li+,Mg+2,Ca+2]")
    if mol.HasSubstructMatch(metal_pattern):
        return False, "Contains metal counter ions - not a pure zwitterion"

    # Look for carboxylate group(s)
    carboxylate_pattern = Chem.MolFromSmarts("[C](=[O])[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not carboxylate_matches:
        return False, "No carboxylate group found"

    # Look for charged nitrogen groups
    nh3_pattern = Chem.MolFromSmarts("[NH3+]")
    nh2_pattern = Chem.MolFromSmarts("[NH2+]")
    nh_pattern = Chem.MolFromSmarts("[NH+]")
    n_pattern = Chem.MolFromSmarts("[N+]")
    
    n_matches = (mol.GetSubstructMatches(nh3_pattern) + 
                mol.GetSubstructMatches(nh2_pattern) +
                mol.GetSubstructMatches(nh_pattern) +
                mol.GetSubstructMatches(n_pattern))
    
    if not n_matches:
        return False, "No positively charged nitrogen group found"

    # Multiple patterns for alpha-amino acid core structure
    core_patterns = [
        # Standard alpha-amino acid
        Chem.MolFromSmarts("[NX4+]-[CH1](-[*])-C(=[O])[O-]"),
        # Cyclic amino acids (like proline)
        Chem.MolFromSmarts("[NX4+]1-[CH2]-[CH2]-[CH1](-C(=[O])[O-])-[CH2]-1"),
        # Alternative connectivity
        Chem.MolFromSmarts("[NX4+]-[CH2]-[CH1](-C(=[O])[O-])-[*]"),
        # Diazaniumyl compounds
        Chem.MolFromSmarts("[NX4+]-[CH1](-[NX3])-C(=[O])[O-]"),
    ]

    found_valid_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(pattern):
            found_valid_core = True
            break

    if not found_valid_core:
        return False, "No valid alpha-amino acid core structure found"

    # Check for reasonable number of charged groups
    n_pos_charges = len(n_matches)
    n_neg_charges = len(carboxylate_matches)
    
    if n_pos_charges > 2 or n_neg_charges > 2:
        return False, "Too many charged groups"
        
    if n_pos_charges != n_neg_charges:
        return False, "Unbalanced charges"

    # Additional validation for dizwitterions
    if n_pos_charges == 2:
        # Check if the charged groups are properly separated
        for match1 in carboxylate_matches:
            for match2 in carboxylate_matches:
                if match1 == match2:
                    continue
                # Get shortest path between carboxylates
                path_len = len(Chem.GetShortestPath(mol, match1[0], match2[0]))
                if path_len < 4:  # Must be reasonably separated
                    return False, "Charged groups too close together"

    return True, "Valid alpha-amino acid zwitterion structure found"