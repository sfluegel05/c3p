"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: CHEBI:36919 thiosugar

A thiosugar is a carbohydrate derivative in which one or more of the oxygens or hydroxy groups 
of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen or any group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sugar backbone
    n_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_rings = rdMolDescriptors.CalcNumRings(mol)

    if n_carbon < 3 or n_oxygen < 3 or n_rings < 1:
        return False, "No valid sugar backbone found"

    # Check for sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
        return False, "No sulfur atoms found"

    # Check if sulfur replaces oxygen or hydroxy group
    for sulfur in sulfur_atoms:
        # Get neighboring atoms
        neighbors = [mol.GetAtomWithIdx(nbr_idx) for nbr_idx in sulfur.GetNeighbors()]

        # Check if sulfur is bound to carbon and oxygen
        if any(nbr.GetAtomicNum() == 6 for nbr in neighbors) and any(nbr.GetAtomicNum() == 8 for nbr in neighbors):
            return True, "Contains sulfur replacing oxygen or hydroxy group in a sugar backbone"

        # Check if sulfur is bound to carbon and sulfur-containing group
        if any(nbr.GetAtomicNum() == 6 for nbr in neighbors) and any(nbr.GetAtomicNum() == 16 for nbr in neighbors):
            return True, "Contains sulfur-containing group replacing oxygen or hydroxy group in a sugar backbone"

        # Check for additional patterns or substructures indicating sulfur replacement (if needed)
        # ...

    return False, "No sulfur found replacing oxygen or hydroxy group in a sugar backbone"