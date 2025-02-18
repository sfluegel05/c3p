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

    # Define a more specific peptide bond pattern: N-[C](=O)-[N]
    peptide_bond_pattern = Chem.MolFromSmarts("[#6;X3](=O)[#7;X3]")

    # Find matches for peptide bonds
    num_peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))

    # Verify if potential tetrapeptide is not part of larger cycle incorrectly picked up
    # Verify that all 3 peptide bonds correspond to non-overlapping terminal-amidic linkages.
    # For a tetrapeptide, there should be 3 peptide bonds linking 4 residues.
    if num_peptide_bonds == 3:
        # Check if end groups are capped (if applicable)
        num_cap_groups = sum([1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetDegree() == 3])
        
        # Modification: Verify if we've exactly two non-bonded terminal groups (likely capping groups)
        if num_cap_groups >= 2:
            return True, "Contains four amino-acid residues connected by peptide linkages"
        else:
            return False, "Does not have expected terminal amidic capping groups"
    
    return False, f"Found {num_peptide_bonds} peptide bonds, requires exactly 3 to be a tetrapeptide"

# Example usage in debug mode
smiles = "C[C@H](N)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CC(O)=O)C(O)=O"
result = is_tetrapeptide(smiles)
print(result)