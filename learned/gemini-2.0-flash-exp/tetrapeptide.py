"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide consists of four amino acid residues connected by three peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the peptide backbone pattern
    peptide_backbone_pattern = Chem.MolFromSmarts("[N]-[C](=[O])-[C]-[N]-[C](=[O])-[C]-[N]-[C](=[O])-[C]-[N]")
    if not mol.HasSubstructMatch(peptide_backbone_pattern):
        return False, "Missing tetrapeptide backbone"

    # 2. Check for peptide bonds (more specific pattern)
    peptide_bond_pattern = Chem.MolFromSmarts("[C](=[O])-[N][C]")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) != 3:
        return False, f"Found {len(peptide_bonds)} peptide bonds, need exactly 3"


    #3. Check for terminal carboxylic acid group. It should be bonded to a carbon.
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX4](=O)O")
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxy_matches) < 1:
        return False, f"Missing terminal carboxylic acid group"
    
    #4. Check for a terminal amine group. It should be bonded to a carbon
    amine_group_pattern = Chem.MolFromSmarts("[CX4][NX3;H2,H1,H0]")
    amine_matches = mol.GetSubstructMatches(amine_group_pattern)
    if len(amine_matches) < 1:
        return False, f"Missing terminal amine group"


    # 5. Check for the correct number of alpha carbons which should be 4.
    alpha_carbon_pattern = Chem.MolFromSmarts("[C;H1]([C])([N])")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    if len(alpha_carbon_matches) != 4:
        return False, f"Found {len(alpha_carbon_matches)} alpha carbons, need exactly 4"

    return True, "Contains four amino acid residues connected by three peptide linkages"