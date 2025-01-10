"""
Classifies: CHEBI:28868 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion typically includes a long aliphatic carbon chain and a deprotonated carboxylate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for deprotonated carboxylate group [O-]=O
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Count the number of carbon atoms to ensure it's a long aliphatic chain (e.g., > 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count <= 8:
        return False, f"Insufficient carbon length for a fatty acid, found {c_count} carbons"

    # Check if the molecular weight falls within typical range for fatty acids
    mol_wt = Descriptors.MolWt(mol)
    if mol_wt < 150:  # Arbitrary lower bound for fatty acid anion (can adjust based on specific cases)
        return False, f"Molecular weight too low, found {mol_wt} Da"

    return True, "Structure matches fatty acid anion with long aliphatic chain and carboxylate group"