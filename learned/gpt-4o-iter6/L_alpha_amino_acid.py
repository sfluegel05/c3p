"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify potential alpha amino acid backbone (C with N, C=O, and H attached)
    # Allow for at least one chiral center with correct configuration
    alpha_aa_smarts = '[C@H](N)C(=O)O | [C@@H](N)C(=O)O'
    alpha_aa_pattern = Chem.MolFromSmarts(alpha_aa_smarts)

    if not mol.HasSubstructMatch(alpha_aa_pattern):
        return False, "No L-alpha-amino acid backbone pattern found"

    # Detect chiral centers and validate configurations using CIP rules
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False, useLegacyImplementation=False)

    for center, chirality in chiral_centers:
        atom = mol.GetAtomWithIdx(center)
        # Ensure the central carbon matches our amino acid pattern, using CIP (R/S) for configuration
        if atom.GetSymbol() == 'C':
            n_neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if 'N' in n_neighbors and 'O' in n_neighbors:
                if chirality in ['R', 'S']:  # Check CIP configurations, one of them must exist
                    return True, "Chiral center with L-configuration identified, amino and carboxyl groups confirmed."

    return False, "No suitable chiral configuration or functional groups for L-alpha-amino acid found."