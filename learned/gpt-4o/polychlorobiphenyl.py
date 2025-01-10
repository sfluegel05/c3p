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
        return False, "Invalid SMILES string"
    
    # Check for biphenyl core (two benzene rings connected by a single bond)
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl core found"
    
    # Count chlorine atoms
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'Cl')

    # Check if chlorine count is between 2 and 10
    if not 2 <= chlorine_count <= 10:
        return False, f"Chlorine atom count is {chlorine_count}, which is not between 2 and 10"

    return True, "Contains a biphenyl core with the correct number of chlorine atoms"

# Examples to test the function
example_smiles = [
    "Clc1ccc(cc1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", # PCB example
    "Clc1ccccc1",  # Not a PCB, single benzene ring
    "Clc1ccccc1-c1ccccc1",  # Biphenyl but only 2 Cl
    "Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)cc(Cl)c1Cl"  # PCB with 6 Cl
]

for smiles in example_smiles:
    is_pcb, reason = is_polychlorobiphenyl(smiles)
    print(f"SMILES: {smiles} -> Is PCB: {is_pcb}, Reason: {reason}")