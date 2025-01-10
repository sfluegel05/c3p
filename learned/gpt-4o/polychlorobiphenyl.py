"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import rdqueries
from rdkit.Chem.rdchem import ChiralType

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
        return False, "Invalid SMILES string"

    # Pattern to match biphenyl core
    biphenyl_pattern = Chem.MolFromSmarts("c1c([a])ccc1-c1c([a])ccc1")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl core found"

    # Count all chlorine atoms
    chlorine_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17])
    if not 2 <= chlorine_count <= 10:
        return False, f"Chlorine count is {chlorine_count}, not within 2-10"

    # Scanning for chlorine positions
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17:  # Check for chlorines
            neighbors = atom.GetNeighbors()
            if not any(neigh.GetIsAromatic() and neigh.GetSymbol() == 'C' for neigh in neighbors):
                return False, "Non-biphenyl or non-aromatic chlorine"

    return True, "Valid polychlorobiphenyl"

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