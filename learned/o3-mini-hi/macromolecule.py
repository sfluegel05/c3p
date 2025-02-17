"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Improved classification for macromolecules.
Definition:
  "A macromolecule is a molecule of high relative molecular mass, 
   the structure of which essentially comprises the multiple repetition 
   of units derived, actually or conceptually, from molecules of low relative molecular mass."

Heuristic improvements:
  1. We count unique occurrences of amide bonds and sugar-like rings
     (using non-overlapping substructure matches).
  2. We apply thresholds that vary with the molecular weight range and total atom count:
         - For weight ≥ 1000 Da: require at least 2 repeating units and at least 40 atoms.
         - For weight between 800 and 1000 Da: require at least 3 repeating units and at least 50 atoms.
         - For weight between 500 and 800 Da: require at least 4 repeating units and at least 60 atoms.
         - Otherwise, classify as not macromolecular.
Usage:
  Call is_macromolecule(smiles) to get a boolean result plus an explanation.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def count_unique_substructures(mol, smarts):
    """
    Count unique, non overlapping occurrences of a substructure in a molecule.
    """
    pattern = Chem.MolFromSmarts(smarts)
    if not pattern:
        return 0
    # Get all substructure matches as sets of atom indices.
    matches = mol.GetSubstructMatches(pattern, uniquify=True)
    # To avoid counting overlapping matches multiple times, we build a set of frozensets.
    unique_matches = set()
    for match in matches:
        unique_matches.add(frozenset(match))
    return len(unique_matches)

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    Improved heuristic:
      - Counts non-overlapping amide bonds ("C(=O)N" pattern) and sugar rings.
      - Applies thresholds that vary with mol weight and total atom count.
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        (bool, str): (True, explanation) if macromolecule criteria are met; otherwise (False, explanation).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute molecular weight and total atom count.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    num_atoms = mol.GetNumAtoms()
    
    # Count unique amide bonds.
    # SMARTS for amide: here we choose a simple pattern; this may catch peptide bonds.
    amide_smarts = "C(=O)N"  
    num_amide = count_unique_substructures(mol, amide_smarts)
    
    # Count unique sugar-like rings.
    # We define a SMARTS that roughly matches pyranoses and can catch furanoses.
    # This is heuristic; you might need to refine as necessary.
    sugar_smarts = "[OX2H][CX4]1[OX2H][CX4][OX2H][CX4][OX2H]1"  
    num_sugar = count_unique_substructures(mol, sugar_smarts)
    
    repetition_evidence = num_amide + num_sugar
    
    # Apply thresholds based on molecular weight.
    if mol_wt >= 1000:
        req_units, req_atoms = 2, 40
    elif 800 <= mol_wt < 1000:
        req_units, req_atoms = 3, 50
    elif 500 <= mol_wt < 800:
        req_units, req_atoms = 4, 60
    else:
        return False, (f"Molecular weight is {mol_wt:.1f} Da with {num_atoms} atoms and "
                       f"{repetition_evidence} repeating units; criteria for a macromolecule are not met.")
    
    if repetition_evidence >= req_units and num_atoms >= req_atoms:
        reason = (f"Molecular weight is {mol_wt:.1f} Da which is high for a macromolecule. "
                  f"Found {num_amide} amide bonds and {num_sugar} sugar units (total evidence: {repetition_evidence}), "
                  f"with {num_atoms} atoms overall. Thresholds: at least {req_units} units and {req_atoms} atoms.")
        return True, reason
    else:
        reason = (f"Molecular weight is {mol_wt:.1f} Da but repetition evidence is low "
                  f"(only {num_amide} amide bonds + {num_sugar} sugar units = {repetition_evidence}) "
                  f"and atom count is {num_atoms}; required thresholds are {req_units} units and {req_atoms} atoms.")
        return False, reason

# Example usage:
if __name__ == '__main__':
    # Example with known macromolecule (peptidic structure)
    test_smiles = "CC[C@H](C)[C@H]1C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)NCCCC[C@@H](C(=O)N[C@@H](C(=O)N1)CCCN)NC(=O)[C@H]([C@@H](C)CC)NC(=O)[C@@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)C2CSC(=N2)C(C(C)CC)N)CC(=O)N)CC(=O)O)CC3=CN=CN3)CC4=CC=CC=C4"
    flag, explanation = is_macromolecule(test_smiles)
    print("Is macromolecule?", flag)
    print("Reason:", explanation)