"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: CHEBI:27815 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a glucoside in which the glycosidic group is derived from D-glucose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for D-glucose substructure
    d_glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)CO")
    d_glucose_matches = mol.GetSubstructMatches(d_glucose_pattern)

    # If no D-glucose substructure found, return False
    if not d_glucose_matches:
        return False, "No D-glucose substructure found"

    # Check if D-glucose is attached to the aglycone via a glycosidic bond
    for match in d_glucose_matches:
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        for atom1, atom2 in zip(atoms, atoms[1:]):
            if atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 8:
                # Found a glycosidic bond (O-C-O)
                break
        else:
            # No glycosidic bond found for this match
            continue
        
        # Check stereochemistry of the glucose moiety
        conf = mol.GetConformer()
        if AllChem.ConfStereochemistry(conf, atoms, mol.GetAtomProperties()):
            # Correct stereochemistry, classify as D-glucoside
            return True, "Contains D-glucose substructure attached to the aglycone via a glycosidic bond"

    # If no match with correct stereochemistry found, return False
    return False, "No D-glucose substructure found with correct stereochemistry and connectivity"