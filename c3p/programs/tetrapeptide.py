"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is any molecule that contains four amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for detecting peptide linkages more accurately
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[CX4][NX3]")  # Includes linkage to next residue
    residue_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)[NX3]") # Updated to include linkage context

    # Match peptide bonds and residues
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    residues = mol.GetSubstructMatches(residue_pattern)

    num_peptide_bonds = len(peptide_bonds)
    num_residues = len(residues)

    # For a tetrapeptide, expect 4 residues typically linked by 3 peptide bonds
    if num_residues == 4 and num_peptide_bonds >= 3:
        return True, "Contains four amino-acid residues connected by peptide linkages"

    return False, f"Found {num_residues} residues and {num_peptide_bonds} peptide bonds, requires 4 residues and 3 peptide bonds"

# Example usage for debug
smiles = "C[C@H](N)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CC(O)=O)C(O)=O"
result = is_tetrapeptide(smiles)
print(result)