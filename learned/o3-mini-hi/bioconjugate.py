"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: A bioconjugate is defined as 'A molecular entity consisting of at least 2 biological molecules covalently linked together'.
Heuristic: We perform a RECAP fragmentation to break the molecule into fragments and then count the fragments 
with at least 6 heavy atoms (as a proxy for biological subunits). If two or more such fragments are found, we classify 
the molecule as a bioconjugate.
"""

from rdkit import Chem
from rdkit.Chem.Recap import RecapDecompose

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    
    A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules 
    covalently linked together. Our heuristic uses RECAP fragmentation to break the molecule into fragments 
    and then counts fragments with at least 6 heavy atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a bioconjugate, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Decompose the molecule using RECAP fragmentation
    try:
        recap_tree = RecapDecompose(mol)
    except Exception as e:
        return False, f"RECAP decomposition failed: {str(e)}"
    
    # Instead of iterating over keys() (which caused the error), use GetAllFragments() to get a dictionary of fragments
    fragments = recap_tree.GetAllFragments()
    if not fragments:
        return False, "Molecule does not decompose into fragments using RECAP; unlikely bioconjugate."
    
    # Count fragments with at least 6 heavy atoms
    bio_fragments_count = 0
    for frag_smiles in fragments.keys():
        frag = Chem.MolFromSmiles(frag_smiles)
        if frag is None:
            continue
        if frag.GetNumHeavyAtoms() >= 6:
            bio_fragments_count += 1

    if bio_fragments_count >= 2:
        return True, f"RECAP fragmentation yielded {bio_fragments_count} fragments with at least 6 heavy atoms."
    else:
        return False, f"Only {bio_fragments_count} biologically significant fragment(s) found."

# Test example
if __name__ == "__main__":
    # Test with gammaGluCys(IAN)Gly
    test_smiles = "N[C@@H](CCC(=O)N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(=O)NCC(O)=O)C(O)=O"
    result, reason = is_bioconjugate(test_smiles)
    print("Is bioconjugate:", result)
    print("Reason:", reason)