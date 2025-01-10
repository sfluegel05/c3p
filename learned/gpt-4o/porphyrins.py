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
    A porphyrin contains a skeleton of four pyrrole nuclei united through the alpha-positions 
    by four methine (-CH=) groups that form a macrocyclic structure, often coordinating a metal ion.

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

    # Define a more flexible SMARTS pattern for porphyrins
    # Pattern description: Four pyrrole rings connected through methine groups, possibly forming a macrocycle
    porphyrin_ring_pattern = Chem.MolFromSmarts("n1c(ccc1)-c2c(-c3[nH]ccc3)-c(-c4[nH]ccc4)n5cccc5")
    
    # Check for the porphyrin macrocycle structure
    if not mol.HasSubstructMatch(porphyrin_ring_pattern):
        return False, "Structure does not match porphyrin macrocycle"

    # Optional: Check presence of metal ions often found in porphyrins
    metal_cations = ["Fe", "Mg", "Zn", "Co", "Ni", "Cu", "Mn"]
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in metal_cations]
    
    if len(metal_atoms) > 0:
        return True, "Contains porphyrin macrocycle with metal coordination"
    else:
        # Return True if the porphyrin structure is detected, even if no metal is present
        return True, "Contains porphyrin macrocycle without metal"

# Example usage:
# is_porphyrins("CC(S)C1=C(C)C2=Cc3c(C(C)S)c(C)c4C=C5C(C)=C(CCC(O)=O)C6=[N+]5[Fe--]5(n34)n3c(=CC1=[N+]25)c(C)c(CCC(O)=O)c3=C6")