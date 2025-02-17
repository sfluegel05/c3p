"""
Classifies: CHEBI:18133 hexose
"""
#!/usr/bin/env python
"""
Classifies: hexose
Definition: Any six‐carbon monosaccharide which in its linear form contains either an aldehyde 
group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).

This function uses two heuristic approaches:
  1. It looks for an open‐chain pattern for aldohexose (or ketohexose) using SMARTS.
  2. If no open‐chain pattern is found it then searches for a cyclic (pyranose‐like) ring: 
     a six‐membered ring with exactly one oxygen (and 5 carbons) whose ring carbons carry several
     hydroxyl (–OH) substituents.
Note that because sugars may appear in various tautomers and as derivatives, this is only an approximate classification.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule qualifies as a hexose.

    A hexose (six‐carbon monosaccharide) is defined here as a molecule that in its open-chain form 
    presents either an aldehyde group at carbon 1 (aldohexose) or a ketone group at carbon 2 (ketohexose), 
    or that contains a cyclic (pyranose) motif with one ring oxygen and several hydroxyl substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a hexose, False otherwise
        str: Explanation for the classification decision
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----- Attempt 1: Look for open-chain hexose patterns -----
    # Aldohexose pattern:
    #   Expected (open-chain) structure: an aldehyde (at carbon 1), followed by 4 secondary alcohol groups, ending with a primary alcohol.
    #   This SMARTS roughly encodes: [CX3H1](=O)[CX4](O)[CX4](O)[CX4](O)[CX4](O)CO
    aldo_hexose_smarts = "[CX3H1](=O)[CX4](O)[CX4](O)[CX4](O)[CX4](O)CO"
    aldo_pattern = Chem.MolFromSmarts(aldo_hexose_smarts)
    if aldo_pattern and mol.HasSubstructMatch(aldo_pattern):
        return True, "Matches open-chain aldohexose pattern (aldehyde at position 1)"
    
    # Ketohexose pattern:
    #   Expected (open-chain) structure for a ketohexose (the ketone at position 2 means carbon chain: HOCH2-C(OH)-C(=O)-...-CH2OH)
    #   This is a heuristic pattern that might capture the core connectivity.
    keto_hexose_smarts = "CO[C;!R](O)C(=O)[CX4](O)[CX4](O)CO"
    keto_pattern = Chem.MolFromSmarts(keto_hexose_smarts)
    if keto_pattern and mol.HasSubstructMatch(keto_pattern):
        return True, "Matches open-chain ketohexose pattern (ketone at position 2)"

    # ----- Attempt 2: Look for a cyclic (pyranose-like) pattern -----
    # Many hexoses exist predominantly in the cyclic (pyranose) form.
    # Here we look for a six-membered ring that contains exactly one oxygen (and 5 carbons).
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            # Count the number of oxygen and carbon atoms in the ring.
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            symbols = [atom.GetSymbol() for atom in ring_atoms]
            if symbols.count("O") == 1 and symbols.count("C") == 5:
                # Look for hydroxyl groups (–OH) attached to the ring carbons.
                hydroxyl_count = 0
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    # Loop over neighbors that are not in the ring.
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() not in ring and nbr.GetSymbol() == "O":
                            # Check if the oxygen is bound to only one heavy atom (typical for an OH group)
                            if nbr.GetDegree() == 1:
                                hydroxyl_count += 1
                # Heuristic: require at least three hydroxyl groups on ring carbons.
                if hydroxyl_count >= 3:
                    return True, "Contains a pyranose ring pattern (six-membered ring with one oxygen and multiple -OH groups)"
    
    return False, "Does not match recognized hexose patterns"

# Example usage:
if __name__ == '__main__':
    # A list of some example SMILES strings (the ones provided are diverse)
    test_smiles = [
        "OC(C(O)CNCCCCCCC)C(O)C(O)CO",  # 1-Deoxy-1-(heptylamino)hexitol
        "O=C(OC1OC(C(O)C(C1O)O)C)C2=CC=CC=C2",  # 1-O-Benzoyl-alpha-L-rhamnopyranoside
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",  # D-allopyranose
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",  # aldehydo-D-galactose (open-chain form)
    ]
    for s in test_smiles:
        result, reason = is_hexose(s)
        print(f"SMILES: {s}\n  Classified as hexose? {result}. Reason: {reason}\n")