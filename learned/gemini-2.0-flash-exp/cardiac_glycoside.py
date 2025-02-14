"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides are characterized by a steroid core, a lactone ring (butenolide), and sugar residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str), True if it is a cardiac glycoside, False otherwise with a reason.
                If parsing fails return (None, None)
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None  # Indicate an issue with SMILES parsing

    # 1. Check for steroid core, specific pattern with ring junctions (more detailed pattern)
    steroid_pattern = Chem.MolFromSmarts("[C]12[C]([C]([C]3)([C]([C]([C]4)([C]1)[C]([C]2)([C]5)[C]([C]3)[C]5)([H]))([H])([H])")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core found."

    # 2. Check for butenolide (5-membered unsaturated lactone) at C17.
    #    Specific attachment at C17: [C]17-[C](=[O])-[C]=[C]-[O]
    butenolide_pattern = Chem.MolFromSmarts("[C]17-[C](=[O])-[C]=[C]-[O]")

    # Find the steroid core match to locate C17
    steroid_match_indices = mol.GetSubstructMatch(steroid_pattern)

    # Check for C17, and find a single carbon which is a neighbour to the steroid core. Then test if the butenolide pattern matches at that carbon
    c17_found = False
    if steroid_match_indices:
        c17_idx_candidates = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and atom.GetIdx() not in steroid_match_indices and all(x.GetIdx() in steroid_match_indices for x in atom.GetNeighbors())]
        
        for c17_idx in c17_idx_candidates:
            # Create SMARTS to check C17 attachment with the butenolide: [C]([C]17)-[C](=[O])-[C]=[C]-[O]
            c17_smarts = f"[C]([{c17_idx}])-[C](=[O])-[C]=[C]-[O]"
            c17_mol = Chem.MolFromSmarts(c17_smarts)

            if c17_mol and mol.HasSubstructMatch(c17_mol):
                 c17_found = True
                 break

    if not c17_found:
      return False, "No butenolide at C17"

    # 3. Check for glycosidic linkage at position C3 of steroid core
    # looking for -[C]-O-[C] where the first carbon is part of the steroid core, and the oxygen is part of glycosidic bond, and the second C is a sugar residue
    glycosidic_pattern = Chem.MolFromSmarts("[C]3-[O]-[C]")
    c3_glycosidic_match_found = False
    if steroid_match_indices:
          c3_idx_candidates = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and all(x.GetIdx() in steroid_match_indices for x in atom.GetNeighbors()) and len([x for x in atom.GetNeighbors() if x.GetAtomicNum() == 8]) > 0 ]
          for c3_idx in c3_idx_candidates:
                c3_smarts = f"[{c3_idx}]-[O]-[C]"
                c3_mol = Chem.MolFromSmarts(c3_smarts)
                if c3_mol and mol.HasSubstructMatch(c3_mol):
                      c3_glycosidic_match_found = True
                      break

    if not c3_glycosidic_match_found:
        return False, "No glycosidic linkage found at C3"
    
    #Additional check for atleast 2 OH groups
    oh_count = 0
    for atom in mol.GetAtoms():
       if atom.GetAtomicNum() == 8 and atom.GetTotalValence() == 2 and len([x for x in atom.GetNeighbors() if x.GetAtomicNum() == 1]) == 1 :
           oh_count += 1
    if oh_count < 2:
       return False, "Too few OH groups. Expected at least two."


    return True, "Meets the criteria for a cardiac glycoside (steroid, butenolide at C17, glycosidic bond at C3, and at least 2 OH groups)."