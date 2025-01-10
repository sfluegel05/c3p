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
    
    # Look specifically for primary ammonium (NH3+)
    nh3_pattern = Chem.MolFromSmarts("[NH3+]")
    nh3_matches = mol.GetSubstructMatches(nh3_pattern)
    
    # Also look for other charged nitrogens for dizwitterions
    other_n_pattern = Chem.MolFromSmarts("[$([NH2+]),$([NH+]),$([N+])]")
    other_n_matches = mol.GetSubstructMatches(other_n_pattern)
    
    if not nh3_matches:
        return False, "No primary ammonium (NH3+) group found"
    
    # Look for carboxylate group (C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[C](=[O])[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not carboxylate_matches:
        return False, "No carboxylate group found"

    # Pattern for alpha-amino acid core structure: NH3+-CH-C(=O)[O-]
    core_pattern = Chem.MolFromSmarts("[NH3+][CH1,CH2][C](=[O])[O-]")
    core_matches = mol.GetSubstructMatches(core_pattern)
    
    # Alternative pattern for diazaniumyl compounds
    diaza_pattern = Chem.MolFromSmarts("[NH3+,NH2+][CH1]([NH1,NH2])[C](=[O])[O-]")
    diaza_matches = mol.GetSubstructMatches(diaza_pattern)
    
    if not (core_matches or diaza_matches):
        return False, "No alpha-amino acid zwitterion core structure found"
    
    # For each potential core, verify it's a proper alpha-amino acid structure
    valid_cores = 0
    for matches, is_diaza in [(core_matches, False), (diaza_matches, True)]:
        for match in matches:
            n_idx = match[0]
            alpha_c_idx = match[1]
            
            # Get the atoms
            alpha_carbon = mol.GetAtomWithIdx(alpha_c_idx)
            n_atom = mol.GetAtomWithIdx(n_idx)
            
            # Check if alpha carbon has correct hybridization (sp3)
            if alpha_carbon.GetHybridization() != Chem.HybridizationType.SP3:
                continue
                
            # Check connectivity
            if alpha_carbon.GetDegree() < 2 or alpha_carbon.GetDegree() > 4:
                continue
                
            # For dizwitterions, allow additional charged groups
            n_charged_groups = len(nh3_matches) + len(other_n_matches)
            n_carboxylates = len(carboxylate_matches)
            
            if n_charged_groups > 2 or n_carboxylates > 2:
                continue
                
            if n_charged_groups != n_carboxylates:
                continue
                
            valid_cores += 1

    if valid_cores > 0:
        return True, "Contains alpha-amino acid zwitterion structure with charged N and COO- groups"
    
    return False, "No valid alpha-amino acid zwitterion structure found"