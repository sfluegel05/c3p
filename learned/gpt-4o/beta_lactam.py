"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam is characterized by a four-membered lactam ring, including an amide nitrogen and a carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure pattern for beta-lactam: a 4-membered ring containing an amide group
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)NC1") # 4-membered ring with amide group
    
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"
    
    return True, "Contains a beta-lactam ring"

# Examples to test if the function works properly (these are beta-lactam structures)
test_smiles = [
    "[H][C@@]1(N)C(=O)N2[C@]1([H])SC(C)(C)[C@]2([H])C([O-])=O",  # 6-aminopenicillanate
    "[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)c1c(C)onc1-c1c(Cl)cccc1Cl)C([O-])=O",  # dicloxacillin(1-)
    "C=1(N2[C@](SCC1/C=C/3\\CCN(C3=O)[C@]4(CN(CC4)C(OCC=5OC(OC5C)=O)=O)[H])([C@@](C2=O)(NC(/C(/C6=NSC(=N6)N)=N\O)=O)[H])[H])C(O)=O",  # ceftobiprole medocaril
]

for smiles in test_smiles:
    result, reason = is_beta_lactam(smiles)
    print(f"SMILES: {smiles} -> Beta-lactam: {result}, Reason: {reason}")