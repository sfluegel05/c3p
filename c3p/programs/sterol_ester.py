"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:37664 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for general sterol backbone (tetracyclic system)
    # More general pattern that matches various sterol configurations
    sterol_pattern = Chem.MolFromSmarts("[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2C)CC=C4[C@@]3(CC[C@@H](C4)[OX2])C")
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol backbone found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester group found"

    # Find oxygen atoms in the sterol core (potential 3-hydroxy position)
    core_oxygens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.IsInRing():
            core_oxygens.append(atom)

    if not core_oxygens:
        return False, "No oxygen in sterol core"

    # Check if any ester oxygen is attached to a core oxygen's carbon
    ester_attached = False
    for ester_match in ester_matches:
        ester_oxygen = mol.GetAtomWithIdx(ester_match[0])
        for core_oxygen in core_oxygens:
            # Check if ester oxygen is connected to the same carbon as core oxygen
            ester_carbon = ester_oxygen.GetNeighbors()[0]
            core_carbon = core_oxygen.GetNeighbors()[0]
            if ester_carbon.GetIdx() == core_carbon.GetIdx():
                ester_attached = True
                break
        if ester_attached:
            break

    if not ester_attached:
        return False, "Ester group not attached to sterol core oxygen"

    # Additional validation: check molecular weight and carbon count
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for sterol ester"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for sterol ester"

    return True, "Contains sterol backbone with ester linkage at core oxygen position"