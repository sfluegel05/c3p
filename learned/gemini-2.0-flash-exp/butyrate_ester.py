"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyryl group (CH3CH2CH2C(=O)-) linked to an oxygen
    that is part of an ester bond (O-C).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the butyrate ester substructure using SMARTS pattern.
    # Matches a butyryl group attached to an ester oxygen linked to any carbon.
    # Also ensures that the carbon is connected to a carbonyl group via an oxygen.
    butyrate_pattern = Chem.MolFromSmarts("CCCC(=O)O[CX3,CX4]")

    #Matches an ester group
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    
    # Check for ester presence. If not, it's not a butyrate ester
    if not mol.HasSubstructMatch(ester_pattern):
      return False, "No ester group found."
    
    # Check to see if a butyrate is connected to the ester
    if mol.HasSubstructMatch(butyrate_pattern):
            # Check to see if it's an anhydride, as this is structurally similar
        anhydride_pattern = Chem.MolFromSmarts("C(=O)O[C](=O)")
        if mol.HasSubstructMatch(anhydride_pattern):
            return False, "Molecule contains an anhydride, not a butyrate ester"

        return True, "Molecule contains a butyrate ester group"
    else:
        return False, "No butyrate ester substructure found"