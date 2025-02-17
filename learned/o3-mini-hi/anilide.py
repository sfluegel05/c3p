"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide – Any aromatic amide obtained by acylation of aniline.
This improved method uses two SMARTS queries to check for a benzene ring directly attached to 
an amide nitrogen. The two queries cover the possibility that the nitrogen is either the substituent
on the benzene (NC(=O)… attached after the ring) or that the benzene ring is attached to the nitrogen 
(NC(=O)… preceded by the ring) as found in many anilide structures.
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is defined as an aromatic amide where the amide nitrogen is derived from aniline –
    that is, the nitrogen is directly attached to a benzene ring (a six-membered aromatic ring made
    solely of carbons) and is also bound to an acyl group.
    
    This improved function uses two SMARTS substructure queries:
      1. "c1ccccc1NC(=O)" – the benzene ring is the substituent through the nitrogen (typical for acylated aniline)
      2. "NC(=O)c1ccccc1" – the acyl group comes first followed by an aromatic ring directly attached to N.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an anilide, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for anilide.
    # Pattern 1: benzene ring connected via an exocyclic amide nitrogen to an acyl group.
    pattern1 = Chem.MolFromSmarts("c1ccccc1NC(=O)")
    # Pattern 2: an amide group whose nitrogen (N) is directly attached to a benzene ring.
    pattern2 = Chem.MolFromSmarts("NC(=O)c1ccccc1")

    if mol.HasSubstructMatch(pattern1):
        return True, "Found anilide substructure (benzene ring attached to N followed by an acyl group)"
    if mol.HasSubstructMatch(pattern2):
        return True, "Found anilide substructure (amide group attached to a benzene ring)"
    
    return False, "No suitable anilide substructure found"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided true positive SMILES.
    test_smiles = "CC(=O)NC1=CC=C(C=C1)C(=S)NCC2=CC=CO2"
    result, reason = is_anilide(test_smiles)
    print("Is anilide:", result)
    print("Reason:", reason)