"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules covalently linked together.
Heuristic: We perform a RECAP fragmentation, which is generally more efficient for large molecules, and then count the number of fragments
with at least 6 heavy atoms (as a proxy for biological subunits). If two or more such fragments exist, we classify the molecule as a bioconjugate.
Note: This heuristic may not capture every nuance.
"""

from rdkit import Chem
from rdkit.Chem.Recap import RecapDecompose
from rdkit.Chem import rdMolDescriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    
    A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules 
    covalently linked together. Our heuristic uses RECAP fragmentation to break the molecule into fragments 
    and counts the fragments having at least 6 heavy atoms.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is likely a bioconjugate, False otherwise
        str: A reason for the classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use RECAP decomposition to fragment the molecule into chemically reasonable parts.
    try:
        recap_tree = RecapDecompose(mol)
    except Exception as e:
        return False, f"RECAP decomposition failed: {str(e)}"

    # If no fragments are returned, then the molecule likely does not consist of separate conjugate parts.
    if not recap_tree:
        return False, "Molecule does not decompose into fragments using RECAP; unlikely bioconjugate."
        
    # The recap_tree is a dictionary where keys are fragment SMILES.
    # Count the number of fragments that are considered 'biological' by having at least 6 heavy atoms.
    bio_fragments_count = 0
    for frag_smiles in recap_tree.keys():
        frag = Chem.MolFromSmiles(frag_smiles)
        if frag is None:
            continue
        if frag.GetNumHeavyAtoms() >= 6:
            bio_fragments_count += 1

    if bio_fragments_count >= 2:
        return True, f"RECAP fragmentation yielded {bio_fragments_count} fragments with at least 6 heavy atoms."
    else:
        return False, f"Only {bio_fragments_count} biologically significant fragment(s) found."
        
# Example usage (try one of the provided SMILES):
if __name__ == "__main__":
    # Example: gammaGluCys(IAN)Gly
    test_smiles = "N[C@@H](CCC(=O)N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(=O)NCC(O)=O)C(O)=O"
    result, reason = is_bioconjugate(test_smiles)
    print("Is bioconjugate:", result)
    print("Reason:", reason)