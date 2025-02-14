"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine has a sphingosine backbone with a fatty acyl group attached to the nitrogen via an amide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for sphingosine backbone, relaxing the original pattern
    # Looking for a chain of carbons with a double bond, two hydroxyl groups, and a nitrogen and CH2OH at the end
    sphingosine_pattern = Chem.MolFromSmarts("[CX4,CX3]-[CX4,CX3](O)-[CX4,CX3](N)-[CX4,CX3](O)-[CX4,CX3]=[CX4,CX3]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Define a SMARTS for a CH2OH group attached to the sphingosine carbon chain
    ch2oh_pattern = Chem.MolFromSmarts("[CX4,CX3]([OX2])")
    ch2oh_matches = mol.GetSubstructMatches(ch2oh_pattern)
    if len(ch2oh_matches) < 1:
        return False, "No CH2OH group found in the sphingosine chain"

    # Define SMARTS for the amide bond
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide bonds, need exactly 1"

    # Define a SMARTS for fatty acid (long carbon chain) attached to the amide bond
    fatty_acyl_pattern = Chem.MolFromSmarts("C[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)

    #Check that the long chain is attached to the amide group
    found_acyl_chain=False
    for match in amide_matches:
        for fatty_match in fatty_acyl_matches:
             amide_nitrogen_atom_index = match[1]
             for atom_index in fatty_match:
                 if mol.GetAtomWithIdx(atom_index).GetBonds()[0].GetOtherAtomIdx(atom_index) == amide_nitrogen_atom_index:
                        found_acyl_chain=True
                        break
             if found_acyl_chain:
                  break

    if not found_acyl_chain:
        return False, "No fatty acyl chain attached to the amide group"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Fatty acid chain too short"

    # Check for at least two hydroxyl groups
    oh_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and len(atom.GetNeighbors()) == 1)
    if oh_count < 2:
        return False, "Must have at least two hydroxyl groups"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
      return False, "Molecular weight too low for N-acylsphingosine"

    return True, "Contains a sphingosine backbone with a fatty acyl group attached via an amide bond"