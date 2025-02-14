"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group and a carboxylic acid group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an alpha-amino acid
    # Handling protonated/deprotonated carboxylic acid and amino groups
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H0,H1,H2][CX4][CX3](=[OX1,OX2])O")
    
    # Check for the substructure match
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
            return False, "No alpha-amino acid substructure found."

    # Additional check to verify that it is a simple alpha amino acid
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    carboxy_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)

    if n_count != 1:
        return False, f"More than one nitrogen found. Expected one for a simple alpha-amino-acid."

    if len(carboxy_matches) != 1:
        return False, f"Expected one carboxy group, found {len(carboxy_matches)}."


    return True, "Molecule matches the alpha-amino acid definition."