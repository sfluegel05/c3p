"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: A bioconjugate is defined as 'A molecular entity consisting of at least 2 biological molecules 
covalently linked together'.

Heuristic: We perform a BRICS fragmentation of the molecule and count the resulting fragments 
with at least 6 heavy atoms as a proxy for biologically meaningful subunits.
"""

from rdkit import Chem
from rdkit.Chem import BRICS

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    
    A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules 
    covalently linked together. Our heuristic uses BRICS fragmentation to break the molecule into fragments
    and then counts fragments with at least 6 heavy atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a bioconjugate, False otherwise.
        str: A reason or explanation for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use BRICS fragmentation which is computationally efficient for complex molecules.
    try:
        fragments = BRICS.BRICSDecompose(mol)
    except Exception as e:
        return False, f"BRICS decomposition failed: {str(e)}"
    
    if not fragments:
        return False, "No fragments generated after BRICS decomposition; unlikely bioconjugate."

    # Count fragments with at least 6 heavy atoms as a proxy for biological subunits.
    bio_fragments_count = 0
    for frag_smiles in fragments:
        frag_mol = Chem.MolFromSmiles(frag_smiles)
        if frag_mol is None:
            continue
        if frag_mol.GetNumHeavyAtoms() >= 6:
            bio_fragments_count += 1

    # We require at least 2 such fragments
    if bio_fragments_count >= 2:
        return True, f"BRICS fragmentation yielded {bio_fragments_count} significant fragments."
    else:
        return False, f"Only {bio_fragments_count} significant fragment(s) were found."

# Example test (can be run when this module is executed directly)
if __name__ == "__main__":
    # For example, testing with gammaGluCys(IAN)Gly
    test_smiles = "N[C@@H](CCC(=O)N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(=O)NCC(O)=O)C(O)=O"
    result, reason = is_bioconjugate(test_smiles)
    print("Is bioconjugate:", result)
    print("Reason:", reason)