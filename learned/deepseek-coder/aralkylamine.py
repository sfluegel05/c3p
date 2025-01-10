"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: CHEBI:50305 aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is an alkylamine in which the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one nitrogen atom (amine group)
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atom found (no amine group)"

    # Check for at least one aromatic ring
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if not aromatic_rings:
        return False, "No aromatic ring found"

    # Check for alkylamine group (nitrogen connected to at least one alkyl group)
    alkylamine_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4]")
    if not mol.HasSubstructMatch(alkylamine_pattern):
        return False, "No alkylamine group found"

    # Check if the alkyl group is substituted by an aromatic group
    # Pattern: Alkyl group connected to an aromatic ring
    aralkyl_pattern = Chem.MolFromSmarts("[CX4][c]")
    if not mol.HasSubstructMatch(aralkyl_pattern):
        return False, "No aromatic group attached to the alkyl group"

    # Verify that the aromatic group is connected to the alkylamine
    # We need to ensure that the aromatic group is part of the same molecule as the alkylamine
    # and that it is connected to the alkyl group of the alkylamine
    for nitrogen in nitrogen_atoms:
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom (alkyl group)
                for carbon_neighbor in neighbor.GetNeighbors():
                    if carbon_neighbor.GetIsAromatic():
                        return True, "Contains an alkylamine group with an aromatic substituent"

    return False, "No aromatic group connected to the alkyl group of the alkylamine"