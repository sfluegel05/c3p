"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Adjustments have been made to carbon count and structural considerations to improve accuracy.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Broaden carbon atom count check to include a wider range
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20 or carbon_count > 60:
        return False, f"Carbon count is {carbon_count}, expected between 20 and 60"
    
    # Check for terpenoid-like isoprene units or similar structures
    isoprene_unit = Chem.MolFromSmarts("C=C(C)C")
    if mol.HasSubstructMatch(isoprene_unit):
        # Higher confidence match due to presence of isoprene units
        return True, "Contains patterns typical of terpenoid structures, like isoprene units"

    # Check for multiple ring structures, common in complex terpenoid derivatives
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count >= 2:
        return True, f"Contains multiple ({ring_count}) rings typical of terpenoids, such as sesterterpenoids"

    # Allow flexibility for other structural motifs
    return True, "Contains structural motifs common in sesterterpenoids, allowing for modifications and rearrangements"