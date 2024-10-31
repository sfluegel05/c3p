from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_amino_acid(smiles: str):
    """
    Determines if a molecule is a beta-amino acid.
    A beta-amino acid has an amino group (-NH2, -NHR, or -NR2) on the beta carbon relative to a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Find carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    # Find amino groups (primary, secondary, or tertiary)
    amino_pattern = Chem.MolFromSmarts('[NX3;H2,H1,H0]')
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    
    if not amino_matches:
        return False, "No amino group found"

    # Pattern for beta-amino acid
    beta_amino_pattern = Chem.MolFromSmarts('[NX3;H2,H1,H0][CH2,CH1,C][CH2,CH1,C]C(=O)[OH]')
    if mol.HasSubstructMatch(beta_amino_pattern):
        # Check if the amino group is exactly at beta position
        matches = mol.GetSubstructMatches(beta_amino_pattern)
        for match in matches:
            amino_n = match[0]
            beta_c = match[1]
            alpha_c = match[2]
            carboxyl_c = match[3]
            
            # Verify the connectivity
            path = Chem.GetShortestPath(mol, amino_n, carboxyl_c)
            if path and len(path) == 4:  # Should be exactly 4 atoms in path for beta position
                amino_atom = mol.GetAtomWithIdx(amino_n)
                h_count = amino_atom.GetTotalNumHs()
                if h_count == 2:
                    amino_type = "primary"
                elif h_count == 1:
                    amino_type = "secondary"
                else:
                    amino_type = "tertiary"
                return True, f"Beta-amino acid found with {amino_type} amino group"

    # Check for cyclic beta-amino acids
    cyclic_beta_pattern = Chem.MolFromSmarts('[NX3;H2,H1,H0;R]~[CH2,CH1,C;R]~[CH2,CH1,C;R]~C(=O)[OH]')
    if mol.HasSubstructMatch(cyclic_beta_pattern):
        return True, "Cyclic beta-amino acid found"

    # Additional check for special cases where the amino group is part of a ring
    # but the carboxylic acid is not
    ring_amino_pattern = Chem.MolFromSmarts('[NX3;H2,H1,H0;R]~[CH2,CH1,C;R]~[CH2,CH1,C]C(=O)[OH]')
    if mol.HasSubstructMatch(ring_amino_pattern):
        return True, "Beta-amino acid with ring-containing amino group found"

    return False, "No amino group found at beta position relative to carboxylic acid"
# Pr=1.0
# Recall=1.0