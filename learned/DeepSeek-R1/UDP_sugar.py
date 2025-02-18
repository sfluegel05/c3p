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
    # Adjusted pattern to match examples: n1ccc(=O)[nH]c(=O)1
    uracil_pattern = Chem.MolFromSmarts("[n]1ccc(=O)[nH]c(=O)1")
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "Missing uracil component"

    # Check for diphosphate linkage (O-P-O-P-O)
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O")
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_matches:
        return False, "No diphosphate linkage found"

    # Verify diphosphate connects uracil component to sugar
    # Find the ribose (connected to uracil) and the sugar via diphosphate
    # Get the two oxygen atoms connecting the diphosphate to other parts
    try:
        diphosphate_atoms = diphosphate_matches[0]
        # The diphosphate connects to ribose (from UDP) and the sugar
        # Check that one side connects to a ribose-like structure (part of UDP)
        # and the other to a sugar (ring with oxygen)
        # Get neighboring atoms of the diphosphate
        phosphate_oxygens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetDegree() == 2]
        connected_to_ribose = False
        connected_to_sugar = False
        for o_idx in phosphate_oxygens:
            neighbor = mol.GetAtomWithIdx(o_idx).GetNeighbors()[0]
            # Check if connected to ribose (part of UDP)
            # Ribose is a 5-membered ring with oxygen, connected to uracil
            ribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](CO)[C@H](O)[C@H]1O")
            if neighbor.GetSymbol() == 'C' and mol.GetAtomWithIdx(neighbor.GetIdx()).GetIsAromatic():
                # Check if connected to uracil's ribose
                ribose_matches = mol.GetSubstructMatches(ribose_pattern)
                if ribose_matches:
                    connected_to_ribose = True
            # Check if connected to sugar (ring with oxygen)
            sugar_ring = Chem.MolFromSmarts("[C@H]1O[C@H](O)CCCC1")  # General sugar ring (5 or 6 membered)
            if neighbor.GetSymbol() == 'C' and mol.GetAtomWithIdx(neighbor.GetIdx()).IsInRing():
                if mol.GetSubstructMatch(sugar_ring):
                    connected_to_sugar = True
        if not (connected_to_ribose and connected_to_sugar):
            return False, "Diphosphate not connecting UDP to sugar"
    except IndexError:
        return False, "Invalid diphosphate connectivity"

    # Check for sugar component (any ring with oxygen, 5 or 6 membered)
    sugar_ring_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]~[C@H]~[C@H]~1")  # Matches 5 or 6 membered rings with oxygen
    if not mol.HasSubstructMatch(sugar_ring_pattern):
        return False, "No sugar ring detected"

    return True, "Contains UDP component connected to sugar via diphosphate linkage"