"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: Tetraterpenoid
Definition: Any terpenoid derived from a tetraterpene. The term includes compounds in which the C40 skeleton 
of the parent tetraterpene has been rearranged or modified by the removal of one or more skeletal atoms.
Heuristic criteria:
  - The molecule is expected to have roughly 30-50 carbon atoms.
  - The molecule is expected to contain a long conjugated polyene chain (here, we search for a substructure of 
    at least three consecutive carbon–carbon double bonds).
Note: This is a heuristic approach and may not correctly classify every borderline case.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    
    A tetraterpenoid is derived from a tetraterpene (nominally a C40 skeleton) that may be rearranged 
    or modified by removal of one or more atoms. In addition, many of these molecules (e.g. carotenoids)
    feature extended conjugated polyene chains.
    
    This function applies two heuristic criteria:
      1. Carbon count is in an expected range (30-50).
      2. Presence of a long conjugated chain – here we search for a polyene substructure with at least
         three consecutive carbon–carbon double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a tetraterpenoid, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms in the molecule
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(carbon_atoms)
    if c_count < 30 or c_count > 50:
        return False, f"Carbon count {c_count} outside expected range (30-50) for tetraterpenoids"
    
    # Look for a long conjugated polyene chain.
    # We use a SMARTS pattern for a chain of three consecutive double bonds:
    # "C=C-C=C-C=C" which represents a polyene fragment.
    polyene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No sufficiently long conjugated polyene chain detected"

    # (Optional) One can also check that the molecular weight is in the expected range (~500 Da for core carotenoid)
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 300 or mw > 800:
        return False, f"Molecular weight {mw:.1f} out of range (300-800 Da) for typical tetraterpenoids"
    
    return True, "Molecule meets carbon count, polyene conjugation, and molecular weight criteria typical for tetraterpenoids"

# Example usage:
if __name__ == "__main__":
    # One of the example SMILES strings provided:
    smiles_example = "CC(C)=CCC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)C=C\\C=C(C)C=C\\C1C(C)=CCCC1(C)C"  # sample SMILES
    result, message = is_tetraterpenoid(smiles_example)
    print("Classified as tetraterpenoid?", result)
    print("Reason:", message)