"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    Sphingomyelin d18:1 is characterized by a sphingosine (d18:1) backbone, a phosphocholine
    head group, and a fatty acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple[bool, str]: True if molecule is a sphingomyelin d18:1, False otherwise, plus the reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Relaxed Sphingosine core pattern (C-C(OH)-C-C(NH)-C)
    sphingosine_core_pattern = Chem.MolFromSmarts("[CX4]-[CX4]([OX2])-[CX4]-[CX4]([NX3])-[CX4]")
    if not mol.HasSubstructMatch(sphingosine_core_pattern):
         return False, "Sphingosine core not found"

    # 2. Phosphocholine head group pattern
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine head group not found"

    # 3. Fatty acyl chain (amide bond) pattern (C(=O)-N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
         return False, "Fatty acyl chain not found"
    
    # 4. Find the long chain of carbons that contains the main core.
    # Find one match for the core
    matches = mol.GetSubstructMatches(sphingosine_core_pattern)
    core_atoms = matches[0]
    
    # Now look for a long chain of carbons. It is required to have 18 atoms and one double bond.
    # We will perform an exploration from a core atom to find the chain.
    
    backbone_atoms = list(core_atoms)
    
    # Find an additional carbon attached to the backbone to initiate the carbon tracing
    found = False
    for atom_idx in core_atoms:
      for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
        if neighbor.GetIdx() not in backbone_atoms and neighbor.GetSymbol() == 'C':
            backbone_atoms.append(neighbor.GetIdx())
            found = True
            break
      if found:
        break

    if not found:
        return False, "Error identifying complete backbone chain."
    
    # Find all the carbons of the backbone
    last_carbon_idx = backbone_atoms[-1]
    while True:
      found = False
      for neighbor in mol.GetAtomWithIdx(last_carbon_idx).GetNeighbors():
        if neighbor.GetIdx() not in backbone_atoms and neighbor.GetSymbol() == 'C':
          backbone_atoms.append(neighbor.GetIdx())
          last_carbon_idx = neighbor.GetIdx()
          found = True
          break
      if not found:
         break

    c_count = sum(1 for atom_idx in backbone_atoms if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
    if c_count != 18:
      return False, f"Incorrect backbone carbon count. Should be 18, got {c_count}"
    
    # Count double bonds within the carbon chain
    double_bond_count = 0
    for bond in mol.GetBonds():
      if bond.GetBondType() == Chem.BondType.DOUBLE:
        if bond.GetBeginAtomIdx() in backbone_atoms and bond.GetEndAtomIdx() in backbone_atoms:
           double_bond_count += 1
    
    if double_bond_count != 1:
        return False, f"Incorrect number of double bonds in backbone, should be 1, got {double_bond_count}"


    return True, "Meets all criteria for sphingomyelin d18:1"