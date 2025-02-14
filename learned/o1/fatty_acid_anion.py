"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: fatty acid anion
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group.

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

    # Check for overall negative charge
    charge = mol.GetFormalCharge()
    if charge >= 0:
        return False, "Molecule does not have a negative charge"

    # Look for carboxylate group (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Get the atoms in the carboxylate group
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    carboxylate_atoms = set()
    for match in carboxylate_matches:
        carboxylate_atoms.update(match)

    # Count total number of carbons excluding carboxylate carbon
    carbon_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetIdx() not in carboxylate_atoms:
            carbon_count +=1

    if carbon_count < 4:
        return False, f"Aliphatic chain too short ({carbon_count} carbons), need at least 4"

    # Check that the chain is aliphatic (non-aromatic)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("a")):
        return False, "Molecule contains aromatic rings"

    return True, "Molecule is a fatty acid anion with a carboxylate group and a long aliphatic chain"