"""
Classifies: CHEBI:75769 B vitamin
"""
from rdkit import Chem

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    B vitamins include Thiamine, Riboflavin, Niacin, Pantothenic acid,
    Pyridoxine, Biotin, Folic acid, and Cobalamin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for key B vitamins
    # Ensure patterns are reviewed for structural accuracy
    patterns = {
        "B1": Chem.MolFromSmarts("CN1C=NC(SC1=CC2=CC=C[NH2+]2)=[NH2+]"),  # Thiamine: Include thiazole and pyrimidine rings
        "B2": Chem.MolFromSmarts("O=C(N=C1N=C(NC(=O)C2=C1N(C=NC2=O)C=CC3=CNC(C=3)=C)[O-])[O-]"),  # Riboflavin: Isoalloxazine rings
        "B3": Chem.MolFromSmarts("c1cc(cnc1)C(=O)O"),  # Nicotinic acid (Niacin): Simplified pyridine structure
        "B5": Chem.MolFromSmarts("CC(C)(CO)C(C(=O)NCCC(=O)O)O"),        # Pantothenic acid
        "B6": Chem.MolFromSmarts("CC1=NC=C(C(C1=O)CO)O"),               # Pyridoxine structure
        "B7": Chem.MolFromSmarts("O=C1N[C@@H]2SCC[C@H]2N1"),             # Biotin: Tetrahydrothieno-structure
        "B9": Chem.MolFromSmarts("Nc1nc2c(cc(C(O)=O)c(=O)n2)n1"),       # Folic acid: Pteridine-based
        "B12": Chem.MolFromSmarts("CNC([Co])C1=NC=A[C@H](A)[N]C1"),      # Cobalamin structure part (focus on large corrin)
    }

    # Check for substructure matches
    for vitamin, pattern in patterns.items():
        if pattern is None:
            return None, "Pattern generation failed."
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches pattern of vitamin {vitamin}"
    
    return False, "Does not match any known B vitamin patterns"