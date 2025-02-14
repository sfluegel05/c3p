"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    Macrolides are characterized by a macrocyclic lactone ring (12+ atoms) derived from polyketides.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find lactone group
    lactone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)

    if not lactone_matches:
          return False, "No lactone group found"
    
    for match in lactone_matches: # check all lactone groups until we find a macrolide

        lactone_atom_indices = [match[0], match[1]] # carbon, oxygen
        
        # Find ring that contains the ester bond
        for ring_info in mol.GetRingInfo().AtomRings():
           if lactone_atom_indices[0] in ring_info and lactone_atom_indices[1] in ring_info:
               ring_size = len(ring_info)
               if ring_size >= 12:
                   return True, f"Macrocyclic lactone ring (size: {ring_size}) found."
                
        #If the loop has completed and none of the rings had >12 atoms, the molecule is not a macrolide.
    return False, "No macrocyclic lactone ring (12+ atoms) found."