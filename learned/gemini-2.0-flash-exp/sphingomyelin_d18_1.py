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

    # 1. Relaxed Sphingosine core pattern (C-C(OH)-C-C(NH)-C=C), including correct stereochemistry
    sphingosine_core_pattern = Chem.MolFromSmarts("[C@H]([OH])[C@H]([N])C=C")
    matches = mol.GetSubstructMatches(sphingosine_core_pattern)
    if not matches:
         return False, "Sphingosine core not found"
    
    # Check for long chain attached to sphingosine core
    core_atoms = matches[0]
    long_chain_found = False
    for core_atom_idx in core_atoms:
       for neighbor in mol.GetAtomWithIdx(core_atom_idx).GetNeighbors():
           if neighbor.GetSymbol() == 'C':
               # Trace the chain from the neighbour, until reaching a dead end. Count carbons
               chain_carbons = [neighbor.GetIdx(),core_atom_idx]
               last_carbon = neighbor.GetIdx()
               while True:
                   found_next = False
                   for next_neighbor in mol.GetAtomWithIdx(last_carbon).GetNeighbors():
                       if next_neighbor.GetIdx() not in chain_carbons and next_neighbor.GetSymbol() == 'C':
                          chain_carbons.append(next_neighbor.GetIdx())
                          last_carbon = next_neighbor.GetIdx()
                          found_next = True
                          break
                   if not found_next:
                      break
               if len(chain_carbons) > 12:
                    long_chain_found = True
                    break

    if not long_chain_found:
          return False, "No long chain attached to sphingosine core"

    # 2. Phosphocholine head group pattern
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine head group not found"

    # 3. Fatty acyl chain (amide bond) pattern (C(=O)-N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
         return False, "Fatty acyl chain not found"

    return True, "Meets all criteria for sphingomyelin d18:1"