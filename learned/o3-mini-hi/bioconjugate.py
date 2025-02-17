"""
Classifies: CHEBI:64985 bioconjugate
"""
"""
Classifies: A bioconjugate is defined as 'A molecular entity consisting of at least 2 biological molecules 
covalently linked together.'

Heuristic: Instead of a full BRICS decomposition (which can be computationally intensive for very large molecules),
we try breaking individual candidate bonds (non-ring, single bonds) one at a time. For each bond break we obtain two fragments.
If at least one cleavage yields two fragments that both contain 6 or more heavy atoms (a rough proxy for a biological subunit),
we consider the molecule a bioconjugate.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    
    A bioconjugate is defined as a molecular entity consisting of at least 2 biological molecules
    covalently linked together. Our heuristic attempts to find a single, cleavable covalent bond (non‐ring, single bond)
    that, upon breakage, produces at least two fragments that are sufficiently large (>= 6 heavy atoms each) to be considered 
    independently "biological" subunits.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is likely a bioconjugate, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Loop over bonds that are candidates for cleavage: single bonds not in rings.
    bond_indices = []
    for bond in mol.GetBonds():
        # Consider only single bonds and bonds that are not part of any ring.
        if bond.GetBondType() == Chem.BondType.SINGLE and not bond.IsInRing():
            bond_indices.append(bond.GetIdx())
    
    # If there are no candidate bonds, then the molecule might be too simple to be a bioconjugate.
    if not bond_indices:
        return False, "No candidate bonds found for cleavage"
    
    # Try each candidate bond
    for bidx in bond_indices:
        try:
            # Fragment the molecule at the bond. The function returns a new mol with dummy atoms where the break occurred.
            frag_mol = Chem.FragmentOnBonds(mol, [bidx])
        except Exception as e:
            continue  # if fragmentation fails for this bond, try the next one
        
        # Get the fragments (as separate molecules) from the fragmented mol.
        frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
        if len(frags) < 2:
            continue
        
        # Count fragments that are "biologically meaningful". Here, we require at least 6 heavy atoms.
        significant_fragments = [f for f in frags if f.GetNumHeavyAtoms() >= 6]
        if len(significant_fragments) >= 2:
            return True, f"Bond index {bidx} yields {len(significant_fragments)} significant fragments (each with ≥6 heavy atoms)."
    
    return False, "No bond cleavage yielded two significant fragments, so the molecule is unlikely to be a bioconjugate."

# For testing purposes, you can run:
if __name__ == "__main__":
    # Example SMILES that was previously used as gammaGluCys(IAN)Gly
    test_smiles = "N[C@@H](CCC(=O)N[C@@H](CSC(C#N)c1c[nH]c2ccccc12)C(=O)NCC(O)=O)C(O)=O"
    result, reason = is_bioconjugate(test_smiles)
    print("Is bioconjugate:", result)
    print("Reason:", reason)