"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    A PCB is a biphenyl compound containing between 2 and 10 chlorine atoms
    attached to the two benzene rings.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorobiphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Verify biphenyl core (two phenyl rings connected by a single bond)
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c1ccccc1")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return (False, "No biphenyl core found")

    # Count all chlorine atoms
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    if not 2 <= chlorine_count <= 10:
        return (False, f"Chlorine atom count is {chlorine_count}, which is not between 2 and 10")

    # Ensure all chlorine atoms are directly connected to biphenyl carbons
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'Cl':
            neighbors = atom.GetNeighbors()
            # Assume chlorines are only on biphenyl aromatic carbons
            if not all(neighbor.GetSymbol() == 'C' and neighbor.GetIsAromatic() for neighbor in neighbors):
                return (False, "Chlorine attached to non-aromatic or non-biphenyl carbon")

    return (True, "Contains a biphenyl core with the correct number of chlorine atoms")

# Example SMILES to test the function
example_smiles = [
    "Clc1ccc(cc1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", # Valid PCB example
    "Clc1ccccc1",  # Not a PCB, just single chlorinated benzene
    "C#Cc1ccc(Cl)cc1", # Alkyne attached, not PCB
    "Clc1cc(Cl)c(cc1)-c1ccc(c(Cl)cc1)Cl"  # Valid PCB with 4 Cl
]

for smiles in example_smiles:
    is_pcb, reason = is_polychlorobiphenyl(smiles)
    print(f"SMILES: {smiles} -> Is PCB: {is_pcb}, Reason: {reason}")