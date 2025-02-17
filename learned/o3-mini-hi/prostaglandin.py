"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: Prostanoid derivatives (prostaglandins)
Definition: Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.
Heuristics used:
1. The molecule must contain a five‐membered ring composed solely of carbon atoms (the prostanoid core).
2. The molecule must contain at least one carboxyl/ester group (the C(=O)O motif) as a sign of the acid functionality.
3. The total number of carbon atoms should be in a range roughly compatible with a C20 backbone (allowing modifications).
If any of these conditions are not met, the function returns False with a reason.
"""

from rdkit import Chem

def is_prostaglandin(smiles: str):
    """
    Determines whether the given SMILES string represents a prostaglandin derivative.
    
    The function uses heuristics:
      - Checks that the molecule contains a five-membered (cyclopentane) carbon ring.
      - Checks that the molecule contains at least one carboxyl/ester functional group
        (a C(=O)O motif).
      - Checks that the number of carbon atoms is in an expected range (roughly C17 to C27).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is likely a prostaglandin derivative, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    # Prostanoic acid based structures generally hover around 20 carbons.
    # We allow some variation (e.g. ester modifications or trim/nor compounds).
    if not (17 <= c_count <= 27):
        return False, f"Carbon count {c_count} outside prostaglandin expected range (17-27)"
    
    # Check for prostanoid cyclopentane core.
    # We look for a five-membered ring and then check if all atoms in that ring are carbons.
    ring_found = False
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) == 5:
            # Check if all atoms in the ring are carbon
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                ring_found = True
                break
    if not ring_found:
        return False, "No 5-membered carbon (cyclopentane) ring found as prostanoid core"
    
    # Check for carboxyl/ester motif.
    # The pattern "C(=O)O" is typical for carboxylic acid (or its esters).
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxyl/ester (C(=O)O) motif found"
    
    return True, "Molecule has a cyclopentane core, carboxyl/ester motif and appropriate carbon count for a prostaglandin derivative"

# For testing purposes uncomment the following:
# test_smiles = [
#     "CCCC[C@H](O)\\C=C\\[C@H]1C=CC(=O)[C@@H]1C\\C=C/CCCC(O)=O",  # prostaglandin A2: should be True
#     "CCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",                    # nonacosanoic acid: should be False
# ]
#
# for s in test_smiles:
#     result, reason = is_prostaglandin(s)
#     print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")