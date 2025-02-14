"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a molecule with a D-glucose moiety linked to another molecule via a glycosidic bond.
    The D-glucose is in its pyranose form.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
         bool: True if the molecule is a D-glucoside, False otherwise.
         str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for D-glucose (pyranose form)
    # Explicitly specifies the stereochemistry at all chiral centers, including the anomeric carbon.
    d_glucose_smarts = "[C@H]1([C@H]([C@@H]([C@H]([C@H](O1)CO)O)O)O)"
    d_glucose_pattern = Chem.MolFromSmarts(d_glucose_smarts)

    #Canonical D-glucose SMILES (for verification of found fragment)
    canonical_d_glucose_smiles = "[C@H]1([C@H]([C@@H]([C@H]([C@@H](O1)CO)O)O)O)"
    
    if d_glucose_pattern is None:
        return None, "Invalid SMARTS pattern for D-glucose"

    # Find all substructure matches of D-glucose
    glucose_matches = mol.GetSubstructMatches(d_glucose_pattern)

    if not glucose_matches:
       return False, "No D-glucose moiety found."

    # Iterate through the glucose matches, verify their configuration, and check for glycosidic bonds.
    for match in glucose_matches:
        
       #Extract the found fragment as a SMILES
        d_glucose_fragment_mol = Chem.MolFragment(mol, match)
        d_glucose_fragment_smiles = Chem.MolToSmiles(d_glucose_fragment_mol)

        #Verify the configuration
        if d_glucose_fragment_smiles != Chem.MolToSmiles(Chem.MolFromSmiles(canonical_d_glucose_smiles)):
            continue # Incorrect stereochemistry

       # Now, look for the glycosidic bond
        for atom_idx in match:
             atom = mol.GetAtomWithIdx(atom_idx)
             
             # Look for Oxygen neighbors.
             for neighbor in atom.GetNeighbors():
                 if neighbor.GetAtomicNum() == 8:
                    is_glycosidic = False
                    
                    # Check if the oxygen is linked to a carbon outside the glucose
                    for neighbor_of_neighbor in neighbor.GetNeighbors():
                          if neighbor_of_neighbor.GetIdx() not in match and neighbor_of_neighbor.GetAtomicNum() == 6 :
                            is_glycosidic = True
                            break
                    if is_glycosidic:
                        return True, "Contains a D-glucose moiety linked via a glycosidic bond."

    return False, "No D-glucose moiety linked via a glycosidic bond found"