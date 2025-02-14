"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:  # Lower bound for triterpenoids + sugars
        return False, "Too few carbons for a triterpenoid saponin"

    # Count oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
      return False, "Too few oxygen atoms for a triterpenoid saponin."

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for a triterpenoid saponin"


    #Check ring count
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 3: #Triterpenoids usually have 3 or more rings
        return False, "Too few rings to be a triterpenoid."

    # Look for triterpenoid core (fused 6-membered rings with some flexibility)
    triterpenoid_core_pattern = Chem.MolFromSmarts("[CX4]1[CX4]2[CX4]3[CX4]4[CX4]1[CX4]5[CX4]2[CX4]3[CX4]45")
    if not mol.HasSubstructMatch(triterpenoid_core_pattern):
        return False, "No triterpenoid core detected"

    # Look for glycosidic bonds (-C-O-C)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[CX4]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_matches) < 1:
         return False, "No glycosidic bond found."
    
    
    return True, "Contains a triterpenoid core with glycosidic bond(s) and sugar(s)."