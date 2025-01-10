"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
from rdkit import Chem

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determine if a molecule is a non-proteinogenic amino acid based on its SMILES string.
    
    A non-proteinogenic amino acid is defined here as one with both amino and carboxyl groups,
    but with non-standard modifications in side chains or additional groups not seen in 
    the 20 standard amino acids.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a non-proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for amino group
    amino_group_pattern = Chem.MolFromSmarts('[NX3H2]')
    # SMARTS pattern for carboxyl group
    carboxyl_group_pattern = Chem.MolFromSmarts('C(=O)O')
    
    # Check for amino group and carboxyl group
    has_amino_group = mol.HasSubstructMatch(amino_group_pattern)
    has_carboxyl_group = mol.HasSubstructMatch(carboxyl_group_pattern)
    
    if not (has_amino_group and has_carboxyl_group):
        return False, "Must contain both amino and carboxyl groups"

    # Count the number of chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral_centers = len(chiral_centers)
    
    # Determine if there are unique side chains
    # A simple heuristic: non-proteinogenic amino acids often have more complex or unusual side chains.
    non_standard_side_chain_criteria = num_chiral_centers > 1 or mol.GetNumAtoms() > 15
    
    if non_standard_side_chain_criteria:
        return True, "Contains unique/modified side chain structures distinguishing it from proteinogenic amino acids"
        
    return False, "Does not have unique/modified side chain features expected in non-proteinogenic amino acids"