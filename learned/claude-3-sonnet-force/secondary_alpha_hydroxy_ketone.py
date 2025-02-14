"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:35654 secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    Secondary alpha-hydroxy ketones have a carbonyl group and a hydroxy group linked
    to the same carbon, which is also connected to an organyl group and a hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find the maximum common substructure with a known secondary alpha-hydroxy ketone
    ref_smiles = "CC(=O)C(O)C"  # Acetoin, a simple secondary alpha-hydroxy ketone
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    mcs = rdFMCS.FindMCS([mol, ref_mol], ringMatchEstimator=lambda x, y: 0)

    # Check if the MCS is a valid secondary alpha-hydroxy ketone substructure
    if mcs.smartsString == "[C&D3](=O)[C&D3](O)[CX4]":
        # Check for additional conditions
        carbonyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 1)
        hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalDegree() == 2)
        
        if carbonyl_count == 1 and hydroxyl_count == 1:
            return True, "Contains a secondary alpha-hydroxy ketone group"
        else:
            return False, "Incorrect number of carbonyl or hydroxyl groups"
    else:
        return False, "No secondary alpha-hydroxy ketone group found"