"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene consists of multiple fused aromatic rings, typically hydrocarbon.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Verify that the molecule contains aromatic rings
    if not mol.GetRingInfo().IsAromatic():
        return False, "Molecule does not contain aromatic rings"
    
    # Compute the number of unique rings and fused rings
    ring_info = mol.GetRingInfo()
    num_aromatic_rings = sum(1 for ring in ring_info.AtomRings() if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))

    # Consider ring fusion by checking if the rings share bonds
    num_fused_rings = 0
    for i, ring1 in enumerate(ring_info.BondRings()):
        for ring2 in ring_info.BondRings():
            if ring1 != ring2 and any(bond in ring1 for bond in ring2):
                num_fused_rings += 1
                break

    # Determine overall polycyclic aromatic system possibility
    if num_aromatic_rings >= 2 and num_fused_rings >= 2:
        return True, "Molecule is a polycyclic arene with multiple fused aromatic rings"
    
    return False, f"Molecule does not meet polycyclic arene criteria with {num_aromatic_rings} aromatic rings and {num_fused_rings} fused"

# Example usage
# print(is_polycyclic_arene("c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c21"))  # Should return True with reason