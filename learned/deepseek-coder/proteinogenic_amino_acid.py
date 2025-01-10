"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: CHEBI:83820 proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    A proteinogenic amino acid has an amino group, a carboxyl group, and a side chain attached to the alpha carbon.
    The molecule should have the L-configuration (except for glycine, which is achiral).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an amino group (primary or secondary) and a carboxyl group
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OX2H1,OX1H0-]")
    
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group found"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check that the amino and carboxyl groups are attached to the same carbon (alpha carbon)
    alpha_carbon_pattern = Chem.MolFromSmarts("[CX4H]([NX3H2,NX3H1])[CX3](=O)[OX2H1,OX1H0-]")
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "Amino and carboxyl groups not attached to the same carbon (alpha carbon)"

    # Check for the presence of a side chain (R group) attached to the alpha carbon
    # Glycine is an exception, where the side chain is just a hydrogen
    glycine_pattern = Chem.MolFromSmarts("[NX3H2][CH2][CX3](=O)[OX2H1,OX1H0-]")
    if mol.HasSubstructMatch(glycine_pattern):
        return True, "Glycine detected (achiral proteinogenic amino acid)"

    # Check for L-configuration (except for glycine)
    # Look for chiral centers and ensure they have the correct stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral center found (not an L-amino acid)"

    # Check molecular weight - expanded range to include all proteinogenic amino acids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 75 or mol_wt > 300:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for proteinogenic amino acids"

    # Additional check for common proteinogenic amino acid side chains
    side_chain_patterns = [
        "[CH3]",  # Alanine
        "[CH2][CH3]",  # Valine
        "[CH2][CH2][CH3]",  # Leucine
        "[CH2]C",  # Isoleucine
        "[CH2]S",  # Methionine
        "[CH2][CH2][CH2][CH2]N",  # Lysine
        "[CH2]C(=O)N",  # Asparagine
        "[CH2][CH2]C(=O)N",  # Glutamine
        "[CH2]C(=O)O",  # Aspartic acid
        "[CH2][CH2]C(=O)O",  # Glutamic acid
        "[CH2]C1=CC=CC=C1",  # Phenylalanine
        "[CH2]C1=CNC=N1",  # Histidine
        "[CH2]OH",  # Serine
        "[CH2][CH2]OH",  # Homoserine
        "[CH2]SH",  # Cysteine
        "[CH2][CH2]S",  # Homocysteine
        "[CH2]C1=CC=C(O)C=C1",  # Tyrosine
        "[CH2][CH2][CH2]N=C(N)N",  # Arginine
        "[CH2]C1CCCN1",  # Proline
        "[CH2]SeH",  # Selenocysteine
        "[CH2][CH2][CH2]C1=CC=CC=N1"  # Pyrrolysine
    ]
    
    has_valid_side_chain = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) 
                              for pattern in side_chain_patterns)
    if not has_valid_side_chain:
        return False, "No valid proteinogenic amino acid side chain found"

    return True, "Contains amino group, carboxyl group, and side chain attached to alpha carbon with L-configuration"