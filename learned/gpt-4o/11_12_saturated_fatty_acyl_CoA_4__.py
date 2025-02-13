"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define comprehensive CoA pattern
    # CoA includes a nucleotide scaffold with base adenine and a diphosphate bridge
    coa_pattern = Chem.MolFromSmarts("NC(=O)[C@H](C(C)(C)COP(O)(=O)OP(OCC1OC(CO)C(O)C1O)n2cnc3ncnc3n2)[C@@H](O)CNC(=O)CCNC(=O)C[CH2]SCC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A group not found"

    # Variable position detection for the saturated 11-12 position in the fatty chain
    chain_length = 0
    acyl_chain = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetDegree() == 1 and mol.GetAtomWithIdx(atom.GetIdx()).GetSymbol() == "C" 
                  and not mol.GetAtomWithIdx(atom.GetIdx()).IsInRing()]
    
    if len(acyl_chain) < 12:
        return False, "Insufficient carbons in acyl chain for 11-12 bond saturation"

    # Checking saturation specifically for the 11,12-resolution based on variable indexing
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in acyl_chain and bond.GetEndAtomIdx() in acyl_chain:
            idx1, idx2 = sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])
            if not bond.IsInRing() and idx2 - idx1 == 1 and mol.GetAtomWithIdx(idx1).GetSymbol() == "C" and mol.GetAtomWithIdx(idx2).GetSymbol() == "C":
                chain_length += 1
                if chain_length == 11 and bond.GetBondType().name != 'SINGLE':
                    return False, "11-12 bond is not saturated"

    return True, "Contains CoA group with a saturated bond from 11-12 in fatty acyl chain"