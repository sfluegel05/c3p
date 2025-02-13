"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:30836 catechin
Members of the class of hydroxyflavan that have a flavan-3-ol skeleton and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for flavan-3-ol skeleton
    flavan_pattern = Chem.MolFromSmarts("[C@H]1[C@@H](c2c(O)cc(O)cc2O[C@H]1c3ccc(O)c(O)c3)O")
    if not mol.HasSubstructMatch(flavan_pattern):
        return False, "No flavan-3-ol skeleton found"

    # Check for benzopyran core
    benzo_pattern = Chem.MolFromSmarts("c1c(O)cc2OCc3cccc(O)c3c2c1")
    if not mol.HasSubstructMatch(benzo_pattern):
        return False, "No benzopyran core found"

    # Check for allowed substituents
    allowed_substituents = ['O', 'C', 'N', 'S', 'P', 'Cl', 'Br', 'F']
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if any(atom not in allowed_substituents for atom in atoms):
        return False, "Contains disallowed substituents"

    # Check for ring system
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:
        return False, "Not enough rings for catechin structure"

    # Check molecular weight - catechins typically <1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, "Molecular weight too high for catechin"

    return True, "Molecule contains flavan-3-ol skeleton and benzopyran core, with allowed substituents and ring system"