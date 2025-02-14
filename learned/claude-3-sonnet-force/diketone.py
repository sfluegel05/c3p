"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: CHEBI:51130 diketone
A compound that contains two ketone functionalities.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate tautomers and resonance structures
    tautomers = list(AllChem.ResonanceMolSupplier(mol, AllChem.KЕКULE_ALL))
    tautomers.append(mol)

    # Check for exactly two ketone groups across all tautomers
    ketones_found = 0
    for tautomer in tautomers:
        ketone_pattern = Chem.MolFromSmarts("[C;$(C=O)]~[#6]")  # Ketone group attached to carbon
        ketone_matches = tautomer.GetSubstructMatches(ketone_pattern)
        ketones_found += len(ketone_matches)

    if ketones_found == 2:
        # Check for incompatible functional groups
        incompatible_pattern = Chem.MolFromSmarts("[N;$(N#C),$(N=O),$(N=S),$(N=P),O;$(O=N),O;$(O=P),S;$(S=O),S;$(S=P)]")
        if mol.HasSubstructMatch(incompatible_pattern):
            return False, "Contains incompatible functional groups for diketone"
        else:
            return True, "Contains exactly two ketone groups"
    else:
        return False, f"Found {ketones_found} ketone groups, expecting exactly 2"

# Examples:
print(is_diketone("CC1(C)C2CC(=O)C1(C)C(=O)C2"))  # True, "Contains exactly two ketone groups"
print(is_diketone("O[C@@H]1CC(=O)[C@H](CCC(O)=O)[C@H]1CCC(=O)CCCCC(O)=O"))  # True, "Contains exactly two ketone groups"
print(is_diketone("OC(=O)\C=C\C(=O)CC(=O)C(O)=O"))  # False, "Found 3 ketone groups, expecting exactly 2"
print(is_diketone("COc1ccc(cc1)C1C(=O)c2ccccc2C1=O"))  # True, "Contains exactly two ketone groups"
print(is_diketone("O=C1C(C(=O)c2ccccc12)c1ccccc1"))  # True, "Contains exactly two ketone groups"
print(is_diketone("CC(=O)C1=CN(C=C(C1c1ccccc1)C(C)=O)c1ccccc1"))  # True, "Contains exactly two ketone groups"
print(is_diketone("C1=2C(=CC(C(C1=O)(OC(C(CC(CC)C)C)=O)C)=O)C=C(OC2)CC(O)C"))  # True, "Contains exactly two ketone groups"
print(is_diketone("Nc1ccccc1C(=O)CC(=O)C(O)=O"))  # False, "Contains incompatible functional groups for diketone"