"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound derived from ammonia by replacing one, two or three
    hydrogen atoms by hydrocarbyl groups

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one nitrogen atom
    has_nitrogen = any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    if not has_nitrogen:
        return False, "No nitrogen atom found"

    # Check for primary, secondary, or tertiary amine (N bonded to H or C)
    primary_amine_pattern = Chem.MolFromSmarts("[N;H2][C]")
    secondary_amine_pattern = Chem.MolFromSmarts("[N;H1]([C])[C]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[N]([C])([C])[C]")
    
    
    has_amine = mol.HasSubstructMatch(primary_amine_pattern) or \
               mol.HasSubstructMatch(secondary_amine_pattern) or \
               mol.HasSubstructMatch(tertiary_amine_pattern)
    
    if not has_amine:
        #Check for cases when N is bonded to H or C. For instance, cyclic amines.
        # This pattern includes case like `NC1CCCCC1` 
        amine_pattern = Chem.MolFromSmarts("[N;H0,H1,H2][!#0]")
        if not mol.HasSubstructMatch(amine_pattern):
            return False, "No amine group found"
        
    # Exclude amides (-C(=O)-N-) using a substructure search to eliminate false positives
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    if mol.HasSubstructMatch(amide_pattern):
       # check that the nitrogen of this amide is not *also* part of an amine.
       for match in mol.GetSubstructMatches(amide_pattern):
         amide_n_index = -1
         for atom_idx in match:
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 7:
                amide_n_index = atom_idx
                break
         if amide_n_index == -1: continue # This should not happen but just in case.
         
         # Check if the amide-N is also an amine, if so, don't exclude.
         n_is_amine = False
         if mol.GetAtomWithIdx(amide_n_index).GetTotalNumHs() == 2: #primary
           n_is_amine = any(mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6 for neighbor in mol.GetAtomWithIdx(amide_n_index).GetNeighbors())
           
         elif mol.GetAtomWithIdx(amide_n_index).GetTotalNumHs() == 1: # secondary
           neighbors = mol.GetAtomWithIdx(amide_n_index).GetNeighbors()
           n_is_amine = sum(1 for neighbor in neighbors if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6) == 2
         
         elif mol.GetAtomWithIdx(amide_n_index).GetTotalNumHs() == 0: # tertiary
          neighbors = mol.GetAtomWithIdx(amide_n_index).GetNeighbors()
          n_is_amine = sum(1 for neighbor in neighbors if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6) == 3
          
         if n_is_amine: continue #If also an amine, continue to the other checks.
         else: return False, "Contains amide group" #Otherwise exclude.

    # Exclude nitro groups (-NO2)
    nitro_pattern = Chem.MolFromSmarts("N(=O)=O")
    if mol.HasSubstructMatch(nitro_pattern):
        return False, "Contains nitro group"
    
    # Exclude nitrile groups (-C#N)
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    if mol.HasSubstructMatch(nitrile_pattern):
        return False, "Contains nitrile group"
        
    # Exclude quaternary ammonium
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+]")
    if mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return False, "Contains quaternary ammonium group"

    # If none of the above exclusion applies, it's an amine.
    return True, "Contains a nitrogen atom bonded to at least one non-hydrogen atom and is not amide nor nitro."