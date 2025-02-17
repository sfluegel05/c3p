"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: Secondary ammonium ion
Definition:
    "An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3."
A secondary ammonium ion arises from a secondary amine (R2NH) upon protonation to form R2NH2+,
so the key feature is a positively charged nitrogen attached to exactly two heavy (non-hydrogen) atoms.
"""

from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    The classification is based on detecting a positively charged nitrogen that is bonded to
    exactly two heavy (non-hydrogen) atoms, which is expected when a secondary amine (R2NH)
    is protonated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains at least one secondary ammonium ion center;
              False otherwise.
        str: Detailed reason for the classification.
    """
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens so that hydrogen counts are correctly represented
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms looking for a nitrogen with a positive formal charge
    for atom in mol.GetAtoms():
        # Check if the atom is nitrogen and has a positive formal charge
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() > 0:
            # Count the number of heavy (non-hydrogen) neighbors
            heavy_neighbor_count = sum(
                1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() != 1
            )
            # For protonated secondary amine (R2NH2+), the nitrogen should have exactly 2 heavy substituents.
            if heavy_neighbor_count == 2:
                return True, ("Found a positively charged nitrogen atom with exactly two "
                              "organic substituents, consistent with a protonated secondary amine.")
    
    return False, "No positively charged nitrogen with two organic substituents was found."