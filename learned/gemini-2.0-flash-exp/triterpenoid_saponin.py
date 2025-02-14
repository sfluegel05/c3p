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

    # Count carbons (Triterpenoids have 30 carbons, + carbons in sugars)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30: # a very low bound
        return False, "Too few carbons for a triterpenoid saponin"

    # Count oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Check molecular weight (typically > 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a triterpenoid saponin"

    # look for the triterpenoid core
    # most triterpenoids have a tetracyclic structure
    # Note that this is quite restrictive
    triterpene_core_pattern = Chem.MolFromSmarts("CC12CCC3(C)C(CCC4(C)C(CCC5C(C)(C)CC[C@H](C)C5)C4)C3C2(C)CC1")
    if not mol.HasSubstructMatch(triterpene_core_pattern):
        # Try a looser pattern without methyls
        triterpene_core_pattern_loose = Chem.MolFromSmarts("C12CCC3C(CCC4C(CCC5CC[C@H](C)C5)C4)C3C2CC1")
        if not mol.HasSubstructMatch(triterpene_core_pattern_loose):
            return False, "No triterpenoid core structure found"

    # Look for glycosidic bonds (-C-O-C, where one C is anomeric carbon of sugar)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-[C@H1](O)")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_matches) < 1: # at least one sugar
         return False, "No glycosidic bond found"

    # Check for hydroxyl groups - triterpenoids and sugars are highly hydroxylated
    oh_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 and len(atom.GetBonds()) == 1)

    if oh_count < 4:
      return False, "Too few hydroxyl groups, unlikely for triterpenoid saponin"

    return True, "Contains triterpenoid core, glycosidic bonds and sufficient hydroxyl groups"