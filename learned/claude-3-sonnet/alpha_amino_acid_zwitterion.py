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

    # Count total charge
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, "Total charge must be 0 for a zwitterion"
    
    # Look for positively charged nitrogen (NH3+, NH2+, N+)
    charged_n_pattern = Chem.MolFromSmarts("[$([NH3+]),$([NH2+]),$([NH+])]")
    charged_n_matches = mol.GetSubstructMatches(charged_n_pattern)
    if not charged_n_matches:
        return False, "No positively charged nitrogen found"
    
    # Look for carboxylate group (C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[C](=[O])[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not carboxylate_matches:
        return False, "No carboxylate group found"

    # Pattern for alpha-amino acid core structure with various N charges:
    # [N+]-C-C(=O)[O-] where C is sp3
    core_pattern = Chem.MolFromSmarts("[$([NH3+]),$([NH2+]),$([NH+]),$([N+])][C;X4][C](=[O])[O-]")
    core_matches = mol.GetSubstructMatches(core_pattern)
    
    if not core_matches:
        return False, "No alpha-amino acid zwitterion core structure found"
    
    # For each potential core, verify it's a proper alpha-amino acid structure
    for match in core_matches:
        n_idx, alpha_c_idx, carboxyl_c_idx = match[0:3]
        
        # Get the alpha carbon atom
        alpha_carbon = mol.GetAtomWithIdx(alpha_c_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check if alpha carbon has correct connectivity
        if alpha_carbon.GetDegree() < 2 or alpha_carbon.GetDegree() > 4:
            continue
            
        # Check if nitrogen has appropriate number of bonds
        if n_atom.GetDegree() > 4:
            continue

        # Count number of carboxylate groups
        carboxylate_count = len(carboxylate_matches)
        if carboxylate_count > 2:  # Allow up to dizwitterions
            return False, "Too many carboxylate groups for simple amino acid zwitterion"
            
        # Count number of charged nitrogens
        charged_n_count = len(charged_n_matches)
        if charged_n_count > 2:  # Allow up to dizwitterions
            return False, "Too many charged nitrogens for simple amino acid zwitterion"

        return True, "Contains alpha-amino acid zwitterion structure with charged N and COO- groups"
    
    return False, "No valid alpha-amino acid zwitterion structure found"