"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine consists of an L-carnitine core with an acyl group attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the L-carnitine core with correct stereochemistry using SMARTS
    # the pattern also allows for both the charged and uncharged forms
    # Note the explicit specification of the C-N bond and the chiral hydrogen

    carnitine_core_smarts = "[C@H]([OX2])(CC(=O)[O-,OH])C[N+](C)(C)C"

    carnitine_core_pattern = Chem.MolFromSmarts(carnitine_core_smarts)

    if not mol.HasSubstructMatch(carnitine_core_pattern):
        return False, "L-carnitine core not found"

    # check for any ester group directly connected to the carnitine core oxygen.

    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")


    matches = mol.GetSubstructMatches(carnitine_core_pattern)
    found_ester = False
    for match in matches:
        chiral_carbon_idx = match[0]
        chiral_carbon = mol.GetAtomWithIdx(chiral_carbon_idx)

        for neighbor in chiral_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                oxygen_idx = neighbor.GetIdx()

                # Check if any ester is attached to the oxygen atom.
                for ester_match in mol.GetSubstructMatches(ester_pattern):
                    if oxygen_idx in ester_match:
                        found_ester = True
                        break
                if found_ester:
                    break

        if found_ester:
           break
    
    if not found_ester:
         return False, "No ester bond found directly linked to the carnitine core"


    # check that the carbon at position 1 is indeed chiral (L configuration)
    matches = mol.GetSubstructMatches(carnitine_core_pattern)
    found_correct_oxygen = False
    for match in matches:
        chiral_carbon_idx = match[0]
        chiral_carbon = mol.GetAtomWithIdx(chiral_carbon_idx)
        
        # now verify that it's indeed an L stereocenter by looking for the hydrogen
        hydrogen_idx = -1
        for neighbor in chiral_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 1:
                hydrogen_idx = neighbor.GetIdx()
                break
        if hydrogen_idx == -1:
            return False, "No hydrogen bound to chiral carbon"

        # Check that the hydrogen is pointing up (CW == L)
        
        if chiral_carbon.GetChiralTag() != Chem.ChiralType.CHI_TETRAHEDRAL_CW:
             return False, "Chiral center not L (or S) configuration"
        
        found_correct_oxygen = True
        break
    
    if not found_correct_oxygen:
        return False, "Chiral carbon not correctly bound to an ester"

    # check the number of nitrogens (sanity check, should be 1 quaternary)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1)
    if n_count != 1:
      return False, f"Number of quaternary nitrogens incorrect, should be 1, found {n_count}"

    return True, "Molecule is an O-acyl-L-carnitine"