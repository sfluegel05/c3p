"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: CHEBI:18085 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is a glucosiduronic acid resulting from the formal condensation
    of any substance with beta-D-glucuronic acid to form a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-D-glucuronic acid substructure pattern
    glucuronic_acid_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)C(=O)O)O)O)O)")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    # Check for the presence of a glycosidic bond (O-C bond between the glucuronic acid and another moiety)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)C(=O)O)O)O)")
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_bond_matches) == 0:
        return False, "No glycosidic bond found"

    # Ensure the molecule contains more than just the glucuronic acid moiety
    # This is a simple check to ensure there is another substructure connected via the glycosidic bond
    # We count the number of atoms not part of the glucuronic acid moiety
    glucuronic_acid_atoms = set()
    for match in mol.GetSubstructMatches(glucuronic_acid_pattern):
        glucuronic_acid_atoms.update(match)
    
    total_atoms = mol.GetNumAtoms()
    if len(glucuronic_acid_atoms) >= total_atoms:
        return False, "Molecule contains only the glucuronic acid moiety"

    # Basic stereochemistry check: Ensure the C1-OH is axial (beta configuration)
    # This is a simplified check and may not catch all cases
    for match in mol.GetSubstructMatches(glucuronic_acid_pattern):
        c1_atom = mol.GetAtomWithIdx(match[0])
        if not c1_atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return False, "Glucuronic acid moiety is not in beta configuration"

    return True, "Contains beta-D-glucuronic acid moiety with a glycosidic bond"

# Example usage:
# smiles = "O[C@H]1[C@@H](O[C@@H]([C@@H](O)[C@@H]1O)C(O)=O)Oc1ccccc1C(O)=O"
# print(is_beta_D_glucosiduronic_acid(smiles))