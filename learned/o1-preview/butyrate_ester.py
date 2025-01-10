"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester

Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is an ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define butyryl group pattern (butyric acid minus the OH)
    butyryl_pattern = Chem.MolFromSmarts("CCCC(=O)")

    # Define ester bond pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)[O,N,S]")

    # Find ester bonds in the molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester groups found"

    # For each ester group, check if the acyl part is butyryl group
    for match in ester_matches:
        carbonyl_c = match[0]  # Carbonyl carbon
        heteroatom = match[2]  # Oxygen, nitrogen, or sulfur

        # Get the acyl fragment (butyryl group)
        acyl_indices = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, carbonyl_c, useHs=True)
        acyl_atoms = set()
        for bond_idx in acyl_indices:
            bond = mol.GetBondWithIdx(bond_idx)
            acyl_atoms.add(bond.GetBeginAtomIdx())
            acyl_atoms.add(bond.GetEndAtomIdx())

        acyl_frag = Chem.PathToSubmol(mol, acyl_indices)

        # Check if the acyl fragment matches the butyryl pattern
        if acyl_frag.HasSubstructMatch(butyryl_pattern):
            return True, "Contains butyrate ester group"

    return False, "Does not contain butyrate ester group"