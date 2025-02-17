"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile acid conjugates.
A bile acid conjugate is defined as any bile acid (a molecule with a steroid nucleus 
and a carboxylic acid side chain) that is conjugated to a functional group which 
confers additional hydrophilicity or charge (such as glycine, taurine/amino acids,
sulfate, glucuronate, sugars, or coenzyme A).
Note: This implementation uses heuristic SMARTS-based filters using rdkit.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    
    This function attempts to confirm two aspects:
    1. The molecule contains a steroid-like (bile acid) core.
       Heuristic: Must contain at least four rings.
    2. The molecule contains a conjugate fragment, as detected by 
       one of several SMARTS patterns (e.g. glycine, taurine, sulfate, glucuronate, or sugar).
       
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a bile acid conjugate, False otherwise.
       str: Reason explaining the classification decision.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Heuristic 1: Check for steroid (bile acid) nucleus.
    # Bile acids are built on a fused ring system of 4 rings.
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 4:
        return False, f"Insufficient rings: found {num_rings}, expected at least 4 for a bile acid core."

    # Heuristic 2: Check for one of the conjugate groups.
    # We will use a list of tuples (name, SMARTS pattern). 
    # (For sugars, here we assume a simple hexose derivative.)
    conjugate_patterns = [
        ("glycine", "[NX3;H2][CH2]C(=O)[O,OH]"),
        ("taurine", "NCCS(=O)(=O)[O,OH]?"),  # sometimes terminal oxygen is not explicit
        ("sulfate", "[OS](=O)(=O)[O,OH]?"),
        ("glucuronate", "OC1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1"),
        ("sugar", "OC1C(O)C(O)C(O)C(O)C1O")
        # Coenzyme A is large and may not be easily defined by a short SMARTS,
        # so it is omitted in this heuristic implementation.
    ]
    
    conjugate_found = False
    conjugate_reason = ""
    for name, smarts in conjugate_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # skip if SMARTS invalid
        if mol.HasSubstructMatch(pattern):
            conjugate_found = True
            conjugate_reason = f"Found conjugate group pattern matching {name}."
            break

    if not conjugate_found:
        return False, "No conjugate fragment detected (expected glycine, taurine, sulfate, glucuronate, or sugar)."
    
    # Additional optional heuristic: A bile acid typically has a carboxylic acid lost by conjugation.
    # For conjugated bile acids, the free carboxyl is turned into an amide or ester.
    # One might wish to check that an amide or ester function is present.
    # For example, look for a C(=O)N (amide) or C(=O)O (ester) group.
    # Here we perform a weak check (if at least one amide is present).
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    ester_pattern   = Chem.MolFromSmarts("C(=O)O")
    if not (mol.HasSubstructMatch(amide_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No amide/ester linkage found that would be expected in a conjugated bile acid."

    return True, f"Contains a steroid nucleus with at least 4 rings and conjugate fragment: {conjugate_reason}"