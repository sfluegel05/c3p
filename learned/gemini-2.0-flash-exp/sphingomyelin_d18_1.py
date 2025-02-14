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

    # 1. Sphingosine core pattern
    sphingosine_core_pattern = Chem.MolFromSmarts("[C]-[C]([OH])-[C]([N]-[C]=O)-[C]=[C]")
    matches = mol.GetSubstructMatches(sphingosine_core_pattern)
    if not matches:
        return False, "Sphingosine core not found"
    
    #Check for long chain on the carbons adjacent to the C=C of the sphingosine core.
    long_chain_found = False
    for match in matches:
      carbon_1 = match[0]
      carbon_5 = match[4]
      if (len(mol.GetAtomWithIdx(carbon_1).GetNeighbors()) > 2) or (len(mol.GetAtomWithIdx(carbon_5).GetNeighbors()) > 2):
            long_chain_found = True
            break
    if not long_chain_found:
         return False, "No long chain attached to sphingosine core"
    

    # 2. Phosphocholine head group pattern
    phosphocholine_pattern = Chem.MolFromSmarts("COP(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine head group not found"

    # 3. Fatty acyl chain (amide bond) pattern
    amide_pattern = Chem.MolFromSmarts("[C]=O[N]-[C]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    if not amide_matches:
        return False, "Fatty acyl chain not found"

    # Verify if amide is attached to the core N
    amide_attached_to_core = False
    for amide_match in amide_matches:
        amide_n_idx = amide_match[1] #amide N index
        for core_match in matches:
           core_n_idx = core_match[2]
           if amide_n_idx == core_n_idx:
              amide_attached_to_core=True
              break
        if amide_attached_to_core:
           break
    
    if not amide_attached_to_core:
         return False, "Fatty acyl chain is not attached to the sphingosine N"
        
    return True, "Meets all criteria for sphingomyelin d18:1"