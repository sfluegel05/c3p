"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules covalently linked together.
Heuristic: We perform a BRICS decomposition (a common fragmentation method) and then count the number of fragments that are 
“biologically meaningful” (here defined as fragments with at least 6 heavy atoms). If two or more such fragments are found, 
we classify the molecule as a bioconjugate.
Note: This is a heuristic approach and may not capture every nuance.
"""

from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    
    A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules 
    covalently linked together. Our heuristic uses BRICS fragmentation to see if the molecule 
    breaks into two or more fragments that are sufficiently large (here, having at least 6 heavy atoms).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is likely a bioconjugate, False otherwise.
        str: A reason for the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use BRICS decomposition to fragment the molecule along chemically sensible bonds.
    try:
        fragments_smiles = BRICS.BRICSDecompose(mol)
    except Exception as e:
        return False, f"BRICS decomposition failed: {str(e)}"

    if not fragments_smiles:
        # If no fragments are returned, then the molecule is not decomposable using BRICS rules.
        return False, "Molecule does not decompose into fragments via BRICS; unlikely bioconjugate."
    
    # Count fragments that are considered to be ‘biological’ using a simple criterion:
    # a fragment with at least 6 heavy (non-hydrogen) atoms.
    bio_fragments_count = 0
    for frag_smiles in fragments_smiles:
        frag = Chem.MolFromSmiles(frag_smiles)
        if frag is None:
            continue
        if frag.GetNumHeavyAtoms() >= 6:
            bio_fragments_count += 1

    if bio_fragments_count >= 2:
        return True, f"BRICS fragmentation yielded {bio_fragments_count} fragments with enough heavy atoms."
    else:
        return False, f"Only {bio_fragments_count} biologically significant fragment(s) found."
        
# Example usage (try one of the provided SMILES):
if __name__ == "__main__":
    test_smiles = "N[C@@H](CCC(=O)N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(=O)NCC(O)=O)C(O)=O"  # gammaGluCys(IAN)Gly
    result, reason = is_bioconjugate(test_smiles)
    print("Is bioconjugate:", result)
    print("Reason:", reason)