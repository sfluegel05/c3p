from rdkit import Chem
from rdkit.Chem import AllChem

def is_C_terminal_proteinogenic_amino_acid_residue(smiles: str):
    """
    Determines if a molecule is a C-terminal proteinogenic amino acid residue.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a C-terminal proteinogenic amino acid residue, False otherwise
        str: Reason for classification
    """
    # Replace * with [*] for RDKit compatibility
    smiles = smiles.replace('*', '[*]')
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Check for presence of amino group and wildcard
    amino_pattern = Chem.MolFromSmarts('[NH1][*]')
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group with wildcard found"

    # Check for alpha carbon between COOH and NH*
    alpha_pattern = Chem.MolFromSmarts('C(=O)(O)[CH1](N[*])')
    if not mol.HasSubstructMatch(alpha_pattern):
        alpha_pattern2 = Chem.MolFromSmarts('O=C(O)[CH1](N[*])')
        if not mol.HasSubstructMatch(alpha_pattern2):
            alpha_pattern3 = Chem.MolFromSmarts('OC(=O)[CH1](N[*])')
            if not mol.HasSubstructMatch(alpha_pattern3):
                return False, "No alpha carbon found connecting COOH and NH*"

    # For proline, check for specific ring structure
    if 'N1' in smiles or '1N' in smiles:
        proline_pattern = Chem.MolFromSmarts('O=C(O)[CH1]1CCCN1[*]')
        if mol.HasSubstructMatch(proline_pattern):
            return True, "Valid proline C-terminal residue"

    # Check chirality - should be S/L configuration (except glycine)
    chiral_centers = Chem.FindMolChiralCenters(mol)
    
    # Glycine has no chiral center
    if len(chiral_centers) == 0:
        glycine_pattern = Chem.MolFromSmarts('C(CN[*])(=O)O')
        if mol.HasSubstructMatch(glycine_pattern):
            return True, "Valid glycine C-terminal residue"
    
    # For other amino acids, check chirality
    if len(chiral_centers) > 0:
        # Get the chiral center connected to both NH and COOH
        alpha_match = mol.GetSubstructMatch(alpha_pattern)
        if not alpha_match:
            alpha_match = mol.GetSubstructMatch(alpha_pattern2)
            if not alpha_match:
                alpha_match = mol.GetSubstructMatch(alpha_pattern3)
        
        alpha_carbon_idx = alpha_match[2] # Index of alpha carbon
        for center in chiral_centers:
            if center[0] == alpha_carbon_idx:
                if '@' not in smiles:
                    return False, "Missing required stereochemistry"
                return True, "Valid C-terminal proteinogenic amino acid residue"

    return True, "Valid C-terminal proteinogenic amino acid residue"
# Pr=1.0
# Recall=0.9090909090909091