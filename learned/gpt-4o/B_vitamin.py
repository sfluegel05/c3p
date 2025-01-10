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

    # SMARTS patterns for some B vitamins:
    patterns = {
        "B1_tetra": Chem.MolFromSmarts("Cc1nccn1Cc1scnc1Cc1ncc(=O)n(CCO)c1"), # Thiamine derivatives
        "B2": Chem.MolFromSmarts("Cc1cc2nc3c(nc(oc=c3)c(C)c2)c(C)c1C"), # Riboflavin
        "B3": Chem.MolFromSmarts("c1cc(C)c(C=O)nc1"), # Nicotinic acid
        "B5": Chem.MolFromSmarts("CCC(C(C(=O)O)O)N"), # Pantothenic acid
        "B6": Chem.MolFromSmarts("CC(O)c1cc(N)c(O)c(C)n1"), # Pyridoxine derivatives
        "B7": Chem.MolFromSmarts("C1=CC2=C(S1)N(F)=V3S(N4C6=CC=C2)CC4(CCC5(C3(CC7)))"), # Biotin (indicative features)
        "B9": Chem.MolFromSmarts("OC(=O)c1ccccc1NCc1ccc(cc1)N"), # Folic acid
        "B12": Chem.MolFromSmarts("c1c[C@@H]2[C@@H]3C@@H]4[C@@H]5CC6CN"), # Cobalamin/related structures
    }

    # Check for substructures
    for vitamin, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Matches pattern of vitamin {vitamin}"
    
    return False, "Does not match any known B vitamin patterns"