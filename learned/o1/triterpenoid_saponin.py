"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
"""

from rdkit import Chem
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

    # Check for triterpenoid core (pentacyclic triterpene)
    # Common triterpenoid skeletons include oleanane, ursane, and lupane
    triterpene_smarts = [
        'C1CCC2(C1)CCC3(C2)CCC4(C3)CCC5(C4)CCCC5',  # oleanane core
        'C1CCC2(C1)CCC3(C2)CCC4(C3)CCCC4',          # ursane core
        'C1CCC2(C1)CCC3(C2)CCC4(C3)CCCC4'           # lupane core
    ]

    has_triterpene_core = False
    for smarts in triterpene_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            has_triterpene_core = True
            break

    if not has_triterpene_core:
        return False, "No triterpenoid core structure found"

    # Check for sugar moieties (pyranose rings)
    sugar_smarts = Chem.MolFromSmarts('[O&R]1[C&R][C&R][C&R][C&R][C&R]1')  # 6-membered ring with oxygen
    sugar_matches = mol.GetSubstructMatches(sugar_smarts)
    if len(sugar_matches) == 0:
        return False, "No sugar moieties found"

    # Check for glycosidic bond between triterpenoid and sugar(s)
    # Look for an ether linkage connecting the triterpenoid core to sugar
    glycosidic_bond_found = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Check for O-C bond where O is connected to sugar ring and C is part of triterpenoid core
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or \
               (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8):
                if (atom1.IsInRingSize(6) and atom2.IsInRing() == False) or \
                   (atom2.IsInRingSize(6) and atom1.IsInRing() == False):
                    glycosidic_bond_found = True
                    break

    if not glycosidic_bond_found:
        return False, "No glycosidic bond connecting triterpenoid and sugar moiety found"

    return True, "Contains triterpenoid core with sugar moiety linked via glycosidic bond"