"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid is an amino acid in which the amino group is located on
    the carbon atom at the position alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxylic acid group pattern: [C](=O)[OH]
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Define amino group pattern: [N;H2,H1;!$(N=*)]
    amino_group_pattern = Chem.MolFromSmarts("[N;!$(N=*)]")  # Exclude imines and other N with double bonds

    # For each carboxylic acid group, check for amino group on alpha carbon
    for carboxyl_match in carboxylic_acid_matches:
        carboxyl_carbon_idx = carboxyl_match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

        # Find alpha carbon (carbon attached to carboxyl carbon)
        alpha_carbons = [atom for atom in carboxyl_carbon.GetNeighbors() if atom.GetAtomicNum() == 6]

        for alpha_carbon in alpha_carbons:
            alpha_carbon_idx = alpha_carbon.GetIdx()

            # Check if alpha carbon has an attached amino group
            attached_to_alpha = [nbr for nbr in alpha_carbon.GetNeighbors()]
            for neighbor in attached_to_alpha:
                if neighbor.GetAtomicNum() == 7:
                    # Check if the nitrogen matches the amino group pattern
                    submol = Chem.PathToSubmol(mol, [neighbor.GetIdx()])
                    if submol.HasSubstructMatch(amino_group_pattern):
                        return True, "Amino group attached to alpha carbon next to carboxyl group"

    return False, "No amino group attached to alpha carbon next to carboxyl group"