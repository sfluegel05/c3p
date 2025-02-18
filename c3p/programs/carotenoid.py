"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: Carotenoid (tetraterpenoid, C40 core typically, with a long conjugated polyene chain)
Note: Retinoids (vitamin A derivatives, typically C20) are excluded.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are tetraterpenoids (derived from psi,psi-carotene, C40 core) that
    feature a long conjugated polyene system. They may be cyclized (carotenes, xanthophylls)
    or modified by oxidation, glycosylation, or similar reactions.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as a carotenoid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string using RDKit. Return False if invalid.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count number of carbon atoms in the molecule.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # For a tetraterpenoid core, one roughly expects around 40 carbons.
    # However, many carotenoids have extra substituents (like sugars) so we accept a broader range.
    if carbon_count < 35:
        return False, f"Too few carbon atoms ({carbon_count}) to be a carotenoid (expected ~40 in the core)"
    
    # Exclude retinoids: They are generally C20 derivatives.
    if carbon_count < 30:
        return False, f"Carbon count ({carbon_count}) is more consistent with retinoids than carotenoids"
    
    # Check for the presence of a long conjugated polyene system.
    # This simple SMARTS pattern requires three alternating double bonds (i.e. four sp2 carbons connected by double bonds with intervening single bonds).
    # Carotenoids usually have a long stretch of conjugation.
    polyene_pattern = Chem.MolFromSmarts("[C;X3]=[C;X3]-[C;X3]=[C;X3]-[C;X3]=[C;X3]")
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No extended conjugated polyene system found"
    
    # Optionally, one could check the total number of rotatable bonds or the molecular weight 
    # but given the large structural variability (e.g., glycosylations, oxidations) we skip these checks.
    
    # If the molecule passes the carbon count and conjugation checks, classify as a carotenoid.
    return True, f"Found {carbon_count} carbon atoms and an extended conjugated polyene system consistent with carotenoids"

# Example usage:
if __name__ == "__main__":
    # Replace the following SMILES string with any carotenoid candidate to test.
    test_smiles = "CC(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C=C/C(=C/C)C)C)C)C)C)C"
    result, reason = is_carotenoid(test_smiles)
    print("Is carotenoid?", result)
    print("Reason:", reason)