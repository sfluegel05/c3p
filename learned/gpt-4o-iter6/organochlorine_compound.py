"""
Classifies: CHEBI:36683 organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound contains at least one carbon-chlorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for C-Cl bond pattern
    organochlorine_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2,CX1]-[Cl]")
    
    if mol.HasSubstructMatch(organochlorine_pattern):
        return True, "Contains at least one carbon-chlorine bond"
    else:
        return False, "No carbon-chlorine bond found"

# Example usage
example_smiles = [
    "[H][C@]12C\\C=C\\C\\C(C)=C\\[C@@]3(C)C=C([C@@H](CC)C[C@]33OC(=O)C(C(=O)[C@@]1(CC)[C@]1([H])CC[C@H](C)[C@]([H])(O[C@H]4C[C@@H](O)[C@H](NC(=O)c5[nH]ccc5Cl)[C@@H](C)O4)[C@@]1([H])C=C2)=C3O)C(O)=O",
    "CCOC(=O)[C@@]1(C)OC(=O)N(C1=O)c1cc(Cl)cc(Cl)c1",
    "CCC(C)n1ncn(-c2ccc(cc2)N2CCN(CC2)c2ccc(OC[C@H]3CO[C@@](Cn4cncn4)(O3)c3ccc(Cl)cc3Cl)cc2)c1=O",
    "Slc1cccc(CC)c1Cl",  # This is not a SMILES but proves detection with an explicit Cl bond
]

for smiles in example_smiles:
    result, reason = is_organochlorine_compound(smiles)
    print(f"SMILES: {smiles} -> Is organochlorine: {result}, Reason: {reason}")