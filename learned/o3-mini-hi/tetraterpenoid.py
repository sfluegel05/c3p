"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: Tetraterpenoid
Definition: Any terpenoid derived from a tetraterpene. The term includes compounds in which the C40 skeleton 
of the parent tetraterpene has been rearranged or modified by the removal of one or more skeletal atoms.
Heuristic criteria (improved):
  - Carbon count: between 20 and 120 (to allow appended groups).
  - Molecular weight: between 200 and 2000 Da.
  - Presence of an extended conjugated polyene chain: we require at least 3 consecutive C=C bonds.
    Instead of a single SMARTS pattern we now try several variants to accommodate differences in stereochemistry.
Note: Being a heuristic, our algorithm may still misclassify borderline cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string using improved heuristic criteria.
    
    Heuristics:
      1. The molecule must contain between 20 and 120 carbon atoms.
      2. The molecular weight must be between 200 and 2000 Da.
      3. The molecule must contain an extended conjugated polyene chain (at least three consecutive C=C bonds).
         To allow for different stereochemical notations we try several SMARTS patterns.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a tetraterpenoid by heuristic, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbon_atoms)
    if c_count < 20 or c_count > 120:
        return False, f"Carbon count {c_count} outside acceptable range (20-120) for tetraterpenoids"
    
    # Check molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 200 or mw > 2000:
        return False, f"Molecular weight {mw:.1f} out of range (200-2000 Da) for typical tetraterpenoids"
    
    # Check for an extended conjugated polyene chain.
    # Instead of using just one SMARTS string, we try several versions to allow different stereochemical notations.
    polyene_smarts_list = [
        "C=C-C=C-C=C",       # plain
        "C=C/C=C/C=C",       # forward slashes
        "C/C=C/C=C/C=C",     # chain with one extra C (allowing for terminal variations)
        "C=C\\C=C\\C=C",     # backslashes only
        "C\\C=C/C=C\\C=C"    # mixed slashes
    ]
    polyene_found = False
    for smarts in polyene_smarts_list:
        polyene_pattern = Chem.MolFromSmarts(smarts)
        if polyene_pattern is None:
            continue
        if mol.HasSubstructMatch(polyene_pattern):
            polyene_found = True
            break
    if not polyene_found:
        return False, "No sufficiently long conjugated polyene chain detected (need at least 3 connected C=C bonds)"
    
    return True, "Molecule meets criteria for tetraterpenoid (carbon count, molecular weight, and polyene chain)"

# Example usage:
if __name__ == "__main__":
    # Try one of the provided example SMILES strings for tetraterpenoids.
    example_smiles = "CO[C@@H]1[C@H](C)O[C@@H](O[C@@H](\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)C=C\\C=C(C)C=C\\C1)C(C)(C)O)[C@@H](OC)[C@@H]1O"
    result, reason = is_tetraterpenoid(example_smiles)
    print("Classified as tetraterpenoid?", result)
    print("Reason:", reason)