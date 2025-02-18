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
    uracil_pattern = Chem.MolFromSmarts("[nH]1ccc(=O)[nH]c(=O)1")
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "Missing uracil component"

    # Check for diphosphate linkage (O-P-O-P-O)
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O")
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_matches:
        return False, "No diphosphate linkage found"

    # Find ribose connected to uracil (part of UDP)
    # Ribose is a 5-membered ring with O, connected to uracil's N1
    uridine_pattern = Chem.MolFromSmarts("[nH]1ccc(=O)[nH]c(=O)1-C2O[C@H](CO[P])[C@H](O)[C@H]2O")
    uridine_matches = mol.GetSubstructMatches(uridine_pattern)
    if not uridine_matches:
        return False, "Uridine (uracil+ribose) not found"

    # Verify diphosphate connects ribose to sugar
    # Get the ribose's CH2O-P group
    ribose_phosphate_pattern = Chem.MolFromSmarts("[C@H]1O[C@H](COP(=O)(O)OP(=O)(O)OC)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(ribose_phosphate_pattern):
        return False, "Ribose not connected to diphosphate"

    # Check sugar component: any ring with oxygen (5 or 6 membered)
    sugar_ring_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]~[C@H]~[C@H]~[C@H]~1")  # 5 or 6 membered
    sugar_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    if not sugar_matches:
        return False, "No sugar ring detected"

    # Ensure the diphosphate connects to the sugar's anomeric oxygen
    # The last O in diphosphate should be connected to the sugar's ring oxygen
    for diphosphate in diphosphate_matches:
        # Get the two P atoms in the diphosphate
        p_atoms = [atom for atom in diphosphate if mol.GetAtomWithIdx(atom).GetSymbol() == 'P']
        if len(p_atoms) != 2:
            continue
        p1, p2 = p_atoms

        # Check P2 (sugar side) is connected to a sugar ring oxygen
        p2_neighbors = mol.GetAtomWithIdx(p2).GetNeighbors()
        for neighbor in p2_neighbors:
            if neighbor.GetSymbol() == 'O':
                o_neighbor = neighbor.GetNeighbors()[0]  # O connected to sugar
                if o_neighbor.GetSymbol() == 'C' and o_neighbor.IsInRing():
                    # Check if this C is part of a sugar ring
                    ring_info = mol.GetRingInfo()
                    for ring in ring_info.AtomRings():
                        if o_neighbor.GetIdx() in ring:
                            # Check if ring has oxygen
                            if any(mol.GetAtomWithIdx(a).GetSymbol() == 'O' for a in ring):
                                return True, "Contains UDP component connected to sugar via diphosphate linkage"

    return False, "Diphosphate not properly connecting UDP to sugar"