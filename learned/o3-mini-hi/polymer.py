"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer is defined as a mixture composed of macromolecules of different kinds,
which may be differentiated by composition, length, degree of branching etc.
This program implements a heuristic for classifying a SMILES string as a polymer.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a given SMILES string represents a polymer—
    i.e. a mixture of distinct macromolecules (with different composition, size, branching, etc.).
    
    Heuristic:
      - Parse the molecule.
      - Split into disconnected fragments.
      - Discard fragments that are too small (likely counterions or salts) using a minimal
        molecular weight (e.g. 150 Da) and requiring at least 5 heavy atoms.
      - Check if at least two remaining fragments exist and that they are chemically distinct.
    
    Args:
        smiles (str): SMILES string for the substance.
        
    Returns:
        bool: True if classified as polymer, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Split molecule into disconnected fragments
    fragments = Chem.GetMolFrags(mol, asMols=True)
    
    # Define minimum criteria for a fragment to be considered a macromolecular component:
    # Here, fragments with molecular weight >= 150 Da and at least 5 heavy atoms.
    def is_large_fragment(frag):
        mw = rdMolDescriptors.CalcExactMolWt(frag)
        heavy_atoms = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)  # exclude hydrogens
        return (mw >= 150) and (heavy_atoms >= 5)
    
    # Filter out small fragments such as simple ions or counterions.
    large_frags = [frag for frag in fragments if is_large_fragment(frag)]
    
    if len(large_frags) < 2:
        # If only one or no “large” fragment is detected,
        # then by our definition, this is not a polymer mixture.
        return False, "Does not meet the criteria for a polymer (fewer than two distinct macromolecular fragments)"
    
    # Compare fragments for distinctiveness using their canonical SMILES.
    frag_smiles = [Chem.MolToSmiles(frag, canonical=True) for frag in large_frags]
    unique_fragments = set(frag_smiles)
    
    if len(unique_fragments) >= 2:
        return True, ("Detected a mixture of at least two distinct macromolecular components "
                      "(each with molecular weight >=150 Da and at least 5 heavy atoms)")
    else:
        return False, ("Multiple large fragments were found but they appear to be chemically identical; "
                       "this is more likely a salt than a polymer mixture")
                       
# Example usage (uncomment to test):
# print(is_polymer("C=1(N(C(N=C(N1)N)(C)C)C=2C=CC(=CC2)Cl)N.C3=CC=CC4=C3C=C(C(=C4CC=5C6=C(C=C(C5O)C(=O)O)C=CC=C6)O)C(=O)O.C=7(N(C(N=C(N7)N)(C)C)C=8C=CC(=CC8)Cl)N"))