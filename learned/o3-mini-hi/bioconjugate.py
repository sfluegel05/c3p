"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: A bioconjugate is defined as 'A molecular entity consisting of at least 2 biological molecules 
covalently linked together'.

Heuristic: We perform a RECAP fragmentation and recursively collect all fragments from the resulting tree.
We then count fragments with at least 6 heavy atoms as a proxy for biological subunits. 
If two or more such fragments are found, we classify the molecule as a bioconjugate.
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
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Attempt RECAP fragmentation
    try:
        recap_tree = RecapDecompose(mol)
    except Exception as e:
        return False, f"RECAP decomposition failed: {str(e)}"
    
    # Helper function to recursively collect all fragment SMILES from the RECAP tree.
    def collect_fragments(node):
        # Each node has an attribute 'smiles' and a dictionary 'children' whose values are Recap nodes.
        fragments = {node.smiles: node}  # include the current node fragment
        # If the node has children then recursively include them
        if hasattr(node, "children") and node.children:
            for child in node.children.values():
                fragments.update(collect_fragments(child))
        return fragments

    fragments_dict = collect_fragments(recap_tree)
    if not fragments_dict:
        return False, "Molecule does not decompose into fragments using RECAP; unlikely bioconjugate."
    
    # Count fragments with at least 6 heavy atoms as a proxy for a biological subunit
    bio_fragments_count = 0
    for frag_smiles in fragments_dict.keys():
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
    # Example: gammaGluCys(IAN)Gly
    test_smiles = "N[C@@H](CCC(=O)N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(=O)NCC(O)=O)C(O)=O"
    result, reason = is_bioconjugate(test_smiles)
    print("Is bioconjugate:", result)
    print("Reason:", reason)