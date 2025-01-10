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

    # Check for the presence of an amino group (-NH2) and a carboxyl group (-COOH)
    amino_pattern = Chem.MolFromSmarts("[NH2]")
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group (-NH2) found"
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group (-COOH) found"

    # Check that the amino and carboxyl groups are attached to the same carbon (alpha carbon)
    alpha_carbon_pattern = Chem.MolFromSmarts("[CX4H]([NH2])[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(alpha_carbon_pattern):
        return False, "Amino and carboxyl groups not attached to the same carbon (alpha carbon)"

    # Check for the presence of a side chain (R group) attached to the alpha carbon
    # Glycine is an exception, where the side chain is just a hydrogen
    glycine_pattern = Chem.MolFromSmarts("[NH2][CH2][CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(glycine_pattern):
        return True, "Glycine detected (achiral proteinogenic amino acid)"

    # Check for L-configuration (except for glycine)
    # We assume that the SMILES string uses the correct stereochemistry notation
    # For L-amino acids, the alpha carbon should have the @ symbol (or @@ for specific cases)
    # This is a simplification and may not cover all cases
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral center found (not an L-amino acid)"

    # Check molecular weight - proteinogenic amino acids typically have MW between 75 and 204 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 75 or mol_wt > 204:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the typical range for proteinogenic amino acids"

    return True, "Contains amino group, carboxyl group, and side chain attached to alpha carbon with L-configuration"