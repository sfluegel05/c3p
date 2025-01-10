"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene (polycyclic aromatic hydrocarbon) consists of multiple aromatic rings fused together.
    
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

    # Count the number of aromatic rings
    ri = mol.GetRingInfo()
    aromatic_ring_count = sum(1 for ring in ri.BondRings() if all(mol.GetBondWithIdx(bond_idx).GetIsAromatic() for bond_idx in ring))
    
    # Check for at least two aromatic rings fused together
    aromatic_fusion = len([set(ring) for ring in ri.BondRings() if all(mol.GetBondWithIdx(bond_idx).GetIsAromatic() for bond_idx in ring)]) >= 2

    if aromatic_ring_count < 2 or not aromatic_fusion:
        return False, "Not enough aromatic rings or insufficient fusion to classify as polycyclic arene"

    # Check for non-carbon, non-hydrogen atoms in aromatic rings
    for ring in ri.BondRings():
        if all(mol.GetBondWithIdx(bond_idx).GetIsAromatic() for bond_idx in ring): # Check if whole ring is aromatic
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() not in [6, 1]:  # Carbon or Hydrogen
                    return False, f"Contains non-carbon, non-hydrogen atoms in aromatic rings: {atom.GetSymbol()}"

    # Ensure planar structure with conjugated pi systems (Hückel's rule)
    num_aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    pi_electron_count = Descriptors.NumAromaticHeterocycles(mol) + 2 * num_aromatic_atoms  # pi electrons from aromatic carbons and heterocycles
    if not pi_electron_count % 4 + 2:
        return False, "Does not follow Hückel's rule (planarity and conjugation check failed)"

    return True, "Molecule is a polycyclic arene with multiple aromatic rings fused together"

# Example usage
# print(is_polycyclic_arene("c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c21"))  # Should return True with reason