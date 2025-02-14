"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the glycerol backbone, phosphocholine and ether/ester groups, while specifying stereochemistry on C2
    # Allow for variation in alkyl/acyl chains.
    pattern = Chem.MolFromSmarts('[C]([OX2]-[C])[C@H]([OX2]-[CX3](=[OX1])-[C])([C]-[OX2])-P(=O)(O)(OCC[N+](C)(C))[O-]')
    match = mol.GetSubstructMatch(pattern)
    
    if match:
       # Check for correct stereochemistry
        glycerol_c2 = match[1]
        glycerol_c2_atom = mol.GetAtomWithIdx(glycerol_c2)
        if not glycerol_c2_atom.HasProp('_Chirality') or glycerol_c2_atom.GetProp('_Chirality') != 'R':
            return False, "Incorrect stereochemistry at glycerol C2"
        
        # Get atom indices for alkyl and acyl chains.
        alkyl_start = match[0]
        acyl_start = match[4]

        # Find the full alkyl and acyl chains, including the connecting oxygen.
        alkyl_mol = Chem.PathToSubmol(mol, Chem.GetShortestPath(mol, alkyl_start, match[1])) #Get the path from start of alkyl to the C2 of glycerol
        acyl_mol  = Chem.PathToSubmol(mol, Chem.GetShortestPath(mol, acyl_start, match[1])) #Get the path from start of acyl to the C2 of glycerol
        

        # Calculate rotatable bonds for alkyl and acyl chains
        alkyl_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(alkyl_mol)
        acyl_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(acyl_mol)
        
        if alkyl_rotatable_bonds < 3 or acyl_rotatable_bonds < 3:
           return False, "Alkyl and/or Acyl chains too short."
        
        return True, "Matches 2-acyl-1-alkyl-sn-glycero-3-phosphocholine criteria"
    else:
      return False, "Does not match 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure"