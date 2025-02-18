"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: CHEBI:15727 UDP-sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar has a UDP nucleotide component attached to a sugar via an anomeric diphosphate linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for uracil component (pyrimidine ring with two oxo groups)
    uracil_pattern = Chem.MolFromSmarts("[nH]1ccc(=O)[nH]c1=O")
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "Missing uracil component"

    # Check for diphosphate linkage (O-P-O-P-O)
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O")
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_matches:
        return False, "No diphosphate linkage found"

    # Verify diphosphate connects two sugar components
    # Find atoms between the two phosphate groups
    try:
        linker = diphosphate_matches[0]
        # Get neighboring atoms of the diphosphate
        phosphate_oxygens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetDegree() == 2]
        sugar_oxygen = None
        for o_idx in phosphate_oxygens:
            neighbor = mol.GetAtomWithIdx(o_idx).GetNeighbors()[0]
            if neighbor.GetSymbol() == 'C' and mol.GetAtomWithIdx(neighbor.GetIdx()).IsInRing():
                sugar_oxygen = o_idx
                break
        if not sugar_oxygen:
            return False, "Diphosphate not connected to a sugar anomeric oxygen"
    except IndexError:
        return False, "Invalid diphosphate connectivity"

    # Check for sugar ring (at least 5 or 6-membered ring with oxygen)
    sugar_ring_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(sugar_ring_pattern):
        return False, "No sugar ring detected"

    return True, "Contains UDP component connected to sugar via diphosphate linkage"