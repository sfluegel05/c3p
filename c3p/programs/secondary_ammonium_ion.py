"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:39143 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has a nitrogen atom with +1 charge bonded to exactly two carbons and two hydrogens (protonated secondary amine).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Count carbon neighbors and total hydrogens (explicit + implicit)
            carbon_count = sum(1 for neighbor in atom.GetNeighbors() 
                             if neighbor.GetAtomicNum() == 6)
            total_h = atom.GetTotalNumHs()
            
            # Check valence (should be 4 for sp3 hybridized N+)
            # Secondary ammonium: 2 carbons + 2 hydrogens = 4 valence
            if (carbon_count == 2 and total_h == 2 and 
                atom.GetTotalValence() == 4 and 
                atom.GetHybridization() == Chem.HybridizationType.SP3):
                return True, "Positively charged nitrogen with two carbon substituents and two hydrogens"
    
    return False, "No secondary ammonium group (N+ with two carbons and two hydrogens) found"