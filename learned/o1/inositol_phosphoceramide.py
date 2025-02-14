"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    An inositol phosphoceramide consists of an inositol residue linked via a phosphodiester bridge
    to a ceramide moiety (sphingoid base linked via an amide bond to a fatty acyl chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify inositol ring (cyclohexane ring with hydroxyl groups on each carbon)
    inositol_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    inositol_match = mol.GetSubstructMatch(inositol_pattern)
    if not inositol_match:
        return False, "No inositol ring found"

    # Identify phosphate group (phosphoric acid ester)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    phosphate_match = mol.GetSubstructMatch(phosphate_pattern)
    if not phosphate_match:
        return False, "No phosphate group found"

    # Identify ceramide moiety (amide connected to long-chain sphingoid base)
    ceramide_pattern = Chem.MolFromSmarts("C(=O)NCCO")
    ceramide_match = mol.GetSubstructMatch(ceramide_pattern)
    if not ceramide_match:
        return False, "No ceramide moiety found"

    # Check connectivity between inositol ring and phosphate group
    inositol_atoms = set(inositol_match)
    phosphorus_idx = phosphate_match[0]  # Index of phosphorus atom in phosphate group
    phosphate_atom = mol.GetAtomWithIdx(phosphorus_idx)

    # Check if phosphate is connected to inositol ring via oxygen
    phosphate_connected_to_inositol = False
    for oxygen in phosphate_atom.GetNeighbors():
        if oxygen.GetAtomicNum() == 8:
            oxygen_idx = oxygen.GetIdx()
            for inositol_atom_idx in inositol_atoms:
                path = Chem.rdmolops.GetShortestPath(mol, oxygen_idx, inositol_atom_idx)
                if len(path) == 2:
                    phosphate_connected_to_inositol = True
                    break
            if phosphate_connected_to_inositol:
                break
    if not phosphate_connected_to_inositol:
        return False, "Phosphate group is not connected to inositol ring"

    # Check if phosphate is connected to ceramide moiety via oxygen
    ceramide_atoms = set(ceramide_match)
    phosphate_connected_to_ceramide = False
    for oxygen in phosphate_atom.GetNeighbors():
        if oxygen.GetAtomicNum() == 8:
            oxygen_idx = oxygen.GetIdx()
            for ceramide_atom_idx in ceramide_atoms:
                path = Chem.rdmolops.GetShortestPath(mol, oxygen_idx, ceramide_atom_idx)
                if path and len(path) <= 5:
                    phosphate_connected_to_ceramide = True
                    break
            if phosphate_connected_to_ceramide:
                break
    if not phosphate_connected_to_ceramide:
        return False, "Phosphate group is not connected to ceramide moiety"

    return True, "Contains inositol phosphoceramide structure"