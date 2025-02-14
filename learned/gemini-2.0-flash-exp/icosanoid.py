"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules arising from the oxidation of C20 EFAs.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for C20 backbone: Look for a chain of at least 18 C atoms, with 20 being more common.
    #   Allow for branching. Must contain at least 1 double bond or ring.
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No C18+ backbone found."

    has_double_bond_or_ring =  mol.HasSubstructMatch(Chem.MolFromSmarts("C=C")) or mol.HasSubstructMatch(Chem.MolFromSmarts("C1CCCC1")) or mol.HasSubstructMatch(Chem.MolFromSmarts("C1CCOC1")) or mol.HasSubstructMatch(Chem.MolFromSmarts("C1CO[C]1"))

    if not has_double_bond_or_ring:
        return False, "No double bonds or rings found in backbone"

    # 2. Check for Carboxylic Acid or Ester Group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    if not (mol.HasSubstructMatch(carboxylic_acid_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "No carboxylic acid or ester group found."
    
    # 3. Check for oxygenation in addition to the ester/acid
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Too few oxygen atoms. Found {o_count}, require at least 2 (for oxidation)"


    # 4. Check for Prostane ring pattern with side chains and oxygenation:
    #   This is a more robust way of identifying the core structure of many icosanoids
    prostane_pattern = Chem.MolFromSmarts("[C]1[C]([C])[C]([O])[C]([C])C1")
    prostane_ring_match = mol.GetSubstructMatches(prostane_pattern)
    if prostane_ring_match:
        # Check for hydroxyl on the cyclopentane ring or other icosanoid structure:
            hydroxyl_on_ring = Chem.MolFromSmarts("[C]1[C]([O])[C][C][C]1")
            if not mol.HasSubstructMatch(hydroxyl_on_ring):
                # Check for THF or epoxide rings for thromboxanes
                thf_ring_pattern = Chem.MolFromSmarts("C1CCOC1")
                epoxide_pattern = Chem.MolFromSmarts("C1OC1")
                if not (mol.HasSubstructMatch(thf_ring_pattern) or mol.HasSubstructMatch(epoxide_pattern)):
                      return False, "Has prostane-like ring but no required oxygenation on the ring"
            return True, "Matches Icosanoid criteria (C20 backbone with prostane like ring with required oxygenation)"

    # 5. Check for Leukotriene-like conjugated system
    # A more specific pattern. Includes a thioether.
    leukotriene_pattern = Chem.MolFromSmarts("[C][S][C][C](O)[C]=C[C]=C[C]=C[C]") # includes a thioether
    conjugated_pattern = Chem.MolFromSmarts("[C]=C[C]=C[C]")  # conjugated double bond system
    if mol.HasSubstructMatch(leukotriene_pattern) or mol.HasSubstructMatch(conjugated_pattern):
      return True, "Matches Icosanoid criteria (C20 backbone with conjugated double bonds and oxygenation)"


    # Check for hydroperoxide
    hydroperoxide_pattern = Chem.MolFromSmarts("[O]-[O]")
    if mol.HasSubstructMatch(hydroperoxide_pattern):
      return True, "Matches Icosanoid criteria (C20 backbone with hydroperoxide)"


    return False, "Does not match Icosanoid criteria"