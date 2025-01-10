"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin contains a skeleton made of four pyrrole rings linked by methine (-CH=) bridges,
    forming a macrocyclic structure potentially coordinating a metal ion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a flexible SMARTS pattern for porphyrins
    # Pattern: Four pyrrole rings connected through methine bridges, arranged in a macrocycle
    porphyrin_ring_pattern = Chem.MolFromSmarts("n1c(ccc1)-c2c(-c3[nH]ccc3)-c(-c4[nH]ccc4)n5ccc(c5*c6[nH]ccc6)")
    
    # Check for the porphyrin macrocyclic structure
    if not mol.HasSubstructMatch(porphyrin_ring_pattern):
        return False, "Structure does not match porphyrin macrocycle"

    # Optional: Check presence of metal ions often coordinated in porphyrins
    metal_cations = ["Fe", "Mg", "Zn", "Co", "Ni", "Cu", "Mn"]
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in metal_cations]
    
    if len(metal_atoms) > 0:
        return True, "Contains porphyrin macrocycle with metal coordination"
    else:
        # Return True if the porphyrin structure is detected, even if no metal is present
        return True, "Contains porphyrin macrocycle without metal"

# Example usage:
# is_porphyrins("CC(S)C1=C(C)C2=Cc3c(C(C)S)c(C)c4C=C5C(C)=C(CCC(O)=O)C6=[N+]5[Fe--]5(n34)n3c(=CC1=[N+]25)c(C)c(CCC(O)=O)c3=C6")